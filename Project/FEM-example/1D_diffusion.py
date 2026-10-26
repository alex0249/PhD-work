import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import ast
import inspect
import textwrap

# ==================== INPUTS ====================

# Mesh and material
N1 = 1000         # Number of elements
R = 1.0           # Particle radius
D = 1.0           # Diffusivity
F = 96485.3329    # Faraday constant

# Initial concentration
c_initial = 1.0

# Current density: positive current adds lithium to the particle.
def i(t):
    return 1.0  # Change this expression to prescribe a time-dependent current.

# Time stepping
dt = 0.01
n_steps = 1000

# Plotting
snapshot_steps = {0, 10, 50, 100, 200, 500, 1000}
figure_size = (7, 4.5)
plot_dpi = 180
plot_filename = "concentration_profiles.png"

# ================= COMMENTED EXAMPLES =================
# Uncomment ONE whole example block to override the inputs above.
# These are illustrative inputs; use consistent units for physical modelling.
# Keep dt = 0.01, R = 1.0 and c_initial = 1.0 for these examples.

# Example 1: constant insertion; compare D = 1.0 against D = 0.01.
# Fast diffusion gives flatter profiles; slow diffusion gives surface accumulation.
# D = 1.0  # Repeat with D = 0.01 and compare the saved plots.
# def i(t):
#     return 0.02 * F
# n_steps = 500
# snapshot_steps = {0, 10, 50, 100, 200, 500}
# plot_filename = "constant_insertion_profiles.png"

# Example 2: insertion pulse followed by rest.
# At t = 1 the pulse ends; the profile then relaxes at constant average.
D = 0.1
def i(t):
    return 10.0 * F if t <= 1.0 else 0.0
n_steps = 1000
snapshot_steps = {0, 0.001, 0.1, 1.0, 2.0, 5.0, 10.0}
# plot_filename = "pulse_rest_profiles.png"

# Example 3: insertion followed by equal-duration extraction.
# The final average returns to its initial value; the profile may remain uneven.
# D = 0.1
# def i(t):
#     return 0.03 * F if t <= 2.0 else -0.03 * F
n_steps = 400
snapshot_steps = {0, 100, 200, 210, 300, 400}
plot_filename = "insertion_extraction_profiles.png"

# ================= END OF INPUTS =================

# Derived mesh quantities
N0 = N1 + 1       # Number of vertices
h = R / N1        # Element length
pi = np.pi
r = np.linspace(0.0, R, N0)
# Each element connects two neighbouring vertices
elements = np.column_stack((
    np.arange(N1),
    np.arange(1, N0),
))
C = np.full((len(r), 1), c_initial)

A = np.zeros((N0, N0))  # Global matrix: N0 × N0
B = np.zeros((N0, 1))  # Global vector: N0 × 1
M = np.zeros((N0, N0))  # Global mass matrix: N0 × N0
def local_mass_matrix(rL, rR):
    h = rR - rL

    M_e = (h / 60) * np.array([
        [12*rL**2 + 6*rL*rR + 2*rR**2,
          3*rL**2 + 4*rL*rR + 3*rR**2],

        [ 3*rL**2 + 4*rL*rR + 3*rR**2,
          2*rL**2 + 6*rL*rR + 12*rR**2]
    ])

    return M_e

for e in range(N1):
    dx_dxi = h/2 #Jacobian
    J_e = h/2 # det
    J_e_inv = 2/h # det inverse
    r_left = r[e]
    r_right = r[e + 1]

    A_e = (
        D * (r_right**3 - r_left**3) / (3 * h**2)
    ) * np.array([[1, -1], [-1, 1]])

    Me = local_mass_matrix(r_left, r_right)

    # Assemble local matrices into global matrices
    A[e,     e]     += A_e[0, 0]
    A[e,     e + 1] += A_e[0, 1]
    A[e + 1, e]     += A_e[1, 0]
    A[e + 1, e + 1] += A_e[1, 1]


    M[e,     e]     += Me[0, 0]
    M[e,     e + 1] += Me[0, 1]
    M[e + 1, e]     += Me[1, 0]
    M[e + 1, e + 1] += Me[1, 1]


# Time stepping
K = M + dt * A

# Store the initial profile; later profiles are saved at snapshot_steps.
profiles = [(0.0, C.copy())]

# Volume average: c_bar = (3/R**3) integral_0^R c(r,t) r**2 dr.
# The mass matrix already integrates the linear FEM field with the r**2 weight.
# Since the shape functions sum to one, column sums of M are the exact weights.
# The common 4*pi factor is divided out of M, A and B in this script.
average_weights = np.asarray(M.sum(axis=0)).ravel() / (R**3 / 3.0)
times = np.arange(n_steps + 1) * dt
average_concentrations = np.empty(n_steps + 1)
current_balance_averages = np.empty(n_steps + 1)
average_concentrations[0] = average_weights @ C[:, 0]
current_balance_averages[0] = average_concentrations[0]

# Time stepping loop
for step in range(n_steps):
    time = (step + 1) * dt
    # Backward Euler evaluates the surface current at the new time.
    B[-1, 0] = R**2 * i(time) / F
    rhs = M @ C + dt * B
    C_new = np.linalg.solve(K, rhs)

    dC_dt = (C_new - C) / dt
    C = C_new

    average_concentrations[step + 1] = average_weights @ C[:, 0]
    # dc_bar/dt = 3*i(t)/(F*R), integrated with the same backward Euler rule.
    current_balance_averages[step + 1] = (
        current_balance_averages[step] + dt * 3.0 * i(time) / (F * R)
    )

    if step + 1 in snapshot_steps:
        profiles.append(((step + 1) * dt, C.copy()))



# Read the active current function so plot labels follow the selected example.
current_tree = ast.parse(textwrap.dedent(inspect.getsource(i)))
current_returns = [node for node in ast.walk(current_tree) if isinstance(node, ast.Return)]
current_expression = "; ".join(ast.unparse(node.value) for node in current_returns)
input_summary = (
    f"R = {R:g} | D = {D:g} | c_initial = {c_initial:g} | F = {F:.8g}\n"
    f"Elements = {N1} | dt = {dt:g} | steps = {n_steps} | end time = {n_steps * dt:g}\n"
    + textwrap.fill(f"i(t) = {current_expression} (positive = insertion)", width=85)
)




def add_input_summary(figure):
    """Reserve a footer below the axes for the active simulation inputs."""
    footer_height = 0.25 + 0.18 * len(input_summary.splitlines())
    figure.set_size_inches(figure_size[0], figure_size[1] + footer_height)
    figure.text(
        0.02, 0.025, input_summary, ha="left", va="bottom", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="whitesmoke", edgecolor="0.8"),
    )
    figure.tight_layout(rect=(0, footer_height / figure.get_figheight(), 1, 1))


# Plot concentration against radius at the saved times.
fig, ax = plt.subplots(figsize=figure_size)
for time, concentration in profiles:
    ax.plot(r, concentration[:, 0], label=f"t = {time:.6g}")

ax.set_xlabel("Radius r")
ax.set_ylabel("Concentration C")
ax.ticklabel_format(axis="y", style="plain", useOffset=False)
ax.set_title("Spherical diffusion: concentration profiles")
ax.grid(True, alpha=0.3)
ax.legend()
add_input_summary(fig)

plot_path = Path(__file__).with_name(plot_filename)
fig.savefig(plot_path, dpi=plot_dpi)
print(f"Saved plot to {plot_path}")
# Plot the volume average at every time step, including a current-balance reference.
average_fig, average_ax = plt.subplots(figsize=figure_size)
average_ax.plot(times, average_concentrations, label="FEM volume average")
average_ax.plot(
    times, current_balance_averages, "--", label="Current balance (backward Euler)"
)
average_ax.set_xlabel("Time t")
average_ax.set_ylabel("Volume-average concentration")
average_ax.ticklabel_format(axis="y", style="plain", useOffset=False)
average_ax.set_title("Spherical diffusion: volume-average concentration")
average_ax.grid(True, alpha=0.3)
average_ax.legend()
add_input_summary(average_fig)

average_plot_path = plot_path.with_name(f"{plot_path.stem}_volume_average.png")
average_fig.savefig(average_plot_path, dpi=plot_dpi)
print(f"Saved volume-average plot to {average_plot_path}")
print(f"Final volume-average concentration: {average_concentrations[-1]:.10g}")
print(
    "Maximum difference from discrete current balance: "
    f"{np.max(np.abs(average_concentrations - current_balance_averages)):.3e}"
)
plt.show()
