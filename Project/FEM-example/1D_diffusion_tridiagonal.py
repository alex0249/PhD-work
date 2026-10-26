import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import ast
import inspect
import textwrap

# ==================== INPUTS ====================

# Units: use SI consistently throughout (m, s, mol, A).
# The values below are illustrative benchmark inputs, not calibrated material data.
# Mesh and material
N1 = 1000         # Number of elements
R = 1.0           # Particle radius [m]
D = 0.1           # Diffusivity [m^2/s], constant for Crank's solution
F = 96485.3329    # Faraday constant [C/mol] = [A*s/mol]

# Initial concentration
c_initial = 0.0   # Initial lithium concentration [mol/m^3]

# Constant inward molar flux: D*dc/dr(R,t) = i(t)/F = j_const.
# This represents a constant flux of lithium ions into the particle.
# Crank section 6.3.5 uses outward flux F0 = -j_const.
# j_const: molar flux [mol/(m^2*s)], positive for insertion.
# i(t): interfacial current density [A/m^2], not total current [A].
# i = F*j assumes one electron per lithium ion.
# For a single fully active sphere, total current I = 4*pi*R^2*i [A].
j_const = 3000.0    # [mol/(m^2*s)]

def i(t):
    """Surface current density [A/m^2]; t is time [s]."""
    return F * j_const

# Time stepping
dt = 0.0001       # Time step [s]
n_steps = 10000

# Plotting
# These are integer STEP NUMBERS, not physical times.
# Physical time = step number * dt. Always include the final step.
# snapshot_steps = {0, 1, 100, 200, 500, 1000, n_steps}
snapshot_steps = {0, 100, 500, 1000, 2000, 5000, 10000}

figure_size = (7, 4.5)
plot_dpi = 180
plot_filename = "constant_flux_comparison.png"

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

# Symmetric tridiagonal matrices: store the diagonal and one off-diagonal.
A_diag = np.zeros(N0)
A_off = np.zeros(N1)
M_diag = np.zeros(N0)
M_off = np.zeros(N1)
B = np.zeros((N0, 1))
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

    # Assemble directly into the nonzero diagonals.
    A_diag[e] += A_e[0, 0]
    A_diag[e + 1] += A_e[1, 1]
    A_off[e] += A_e[0, 1]
    M_diag[e] += Me[0, 0]
    M_diag[e + 1] += Me[1, 1]
    M_off[e] += Me[0, 1]


def tridiagonal_matvec(diagonal, off_diagonal, values):
    """Multiply a symmetric tridiagonal matrix by a column vector."""
    result = diagonal[:, None] * values
    result[:-1] += off_diagonal[:, None] * values[1:]
    result[1:] += off_diagonal[:, None] * values[:-1]
    return result


def factor_tridiagonal(diagonal, off_diagonal):
    """Compute K = L D L.T for symmetric positive-definite K.

    D_diag stores the diagonal of D; L_subdiag stores the subdiagonal of L.
    L has ones on its diagonal. Both factors are stored as 1D arrays.
    """
    D_diag = diagonal.copy()
    L_subdiag = np.empty_like(off_diagonal)
    for row in range(len(off_diagonal)):
        L_subdiag[row] = off_diagonal[row] / D_diag[row]
        D_diag[row + 1] -= L_subdiag[row] * off_diagonal[row]
    if not np.all(np.isfinite(D_diag)) or np.any(D_diag <= 0):
        raise ValueError("The diffusion matrix must be positive definite.")
    return D_diag, L_subdiag


def solve_tridiagonal(D_diag, L_subdiag, rhs):
    """Solve K x = rhs using D_diag and L_subdiag from K = L D L.T."""
    solution = rhs[:, 0].copy()
    # Forward substitution: L y = rhs.
    for row in range(1, len(D_diag)):
        solution[row] -= L_subdiag[row - 1] * solution[row - 1]
    # Diagonal solve: D z = y.
    solution /= D_diag
    # Back substitution: L.T x = z.
    for row in range(len(D_diag) - 2, -1, -1):
        solution[row] -= L_subdiag[row] * solution[row + 1]
    return solution[:, None]

def calculate_crank_roots(number_of_terms):
    """Find the positive roots beta_n of tan(beta_n) = beta_n."""

    roots = np.zeros(number_of_terms)

    for n in range(1, number_of_terms + 1):

        # The nth positive root lies in this interval.
        left = n * np.pi
        right = n * np.pi + np.pi / 2

        # Repeatedly halve the interval containing the root.
        for iteration in range(60):

            middle = (left + right) / 2

            value_left = np.sin(left) - left * np.cos(left)
            value_middle = np.sin(middle) - middle * np.cos(middle)

            if value_left * value_middle > 0:
                left = middle
            else:

                right = middle

        roots[n - 1] = (left + right) / 2

    return roots
number_of_terms = 80
# Very early times require more modes to cancel the spatial polynomial accurately.
# Use enough terms that omitted modes have decay below about exp(-36) at t=dt.
# With R=1, D=0.1 and dt=0.0001 this increases 80 terms to 604.
number_of_terms = max(
    number_of_terms, int(np.ceil(6 * R / (np.pi * np.sqrt(D * dt))))
)
beta_roots = calculate_crank_roots(number_of_terms)
print(f"Crank reference terms: {number_of_terms}")

def crank_constant_flux(r, time, R, D, c0, j, beta_roots):
    """Evaluate Crank's constant-flux solution at the supplied radii."""

    r = np.asarray(r, dtype=float)
    x = r / R

    # The initial condition is known exactly.
    if time == 0:
        return np.full_like(r, c0)

    # Accumulate the transient series, one mode at a time.
    series = np.zeros_like(r)

    for beta in beta_roots:

        # g_n(x) = sin(beta*x) / (x*sin(beta)).
        # Using sinc also handles the centre x = 0 safely.
        spatial_mode = (
            np.sinc(beta * x / np.pi)
            / np.sinc(beta / np.pi)
        )

        decay = np.exp(-D * beta**2 * time / R**2)

        series += spatial_mode * decay / beta**2

    # The volume-average concentration increases with insertion.
    mean_concentration = c0 + 3 * j * time / R

    # Spatial variation relative to the volume average.
    spatial_correction = (j * R / D) * (
        x**2 / 2 - 3 / 10 - 2 * series
    )

    return mean_concentration + spatial_correction

# K is constant: factor it once, then reuse its factors at every time step.
K_diag = M_diag + dt * A_diag
K_off = M_off + dt * A_off
K_D_diag, K_L_subdiag = factor_tridiagonal(K_diag, K_off)

# Store the initial profile; later profiles are saved at snapshot_steps.
profiles = [(0.0, C.copy())]

# Volume average: c_bar = (3/R**3) integral_0^R c(r,t) r**2 dr.
# The mass matrix already integrates the linear FEM field with the r**2 weight.
# Since the shape functions sum to one, column sums of M are the exact weights.
# The common 4*pi factor is divided out of M, A and B in this script.
mass_column_sums = M_diag.copy()
mass_column_sums[:-1] += M_off
mass_column_sums[1:] += M_off
average_weights = mass_column_sums / (R**3 / 3.0)
times = np.arange(n_steps + 1) * dt
average_concentrations = np.empty(n_steps + 1)
current_balance_averages = np.empty(n_steps + 1)
# Record the outermost FEM node at EVERY step, not just the profile snapshots.
surface_concentrations = np.empty(n_steps + 1)
surface_concentrations[0] = C[-1, 0]
average_concentrations[0] = average_weights @ C[:, 0]
current_balance_averages[0] = average_concentrations[0]

# Time stepping loop
for step in range(n_steps):
    time = (step + 1) * dt
    # Backward Euler evaluates the surface current at the new time.
    B[-1, 0] = R**2 * i(time) / F
    rhs = tridiagonal_matvec(M_diag, M_off, C) + dt * B
    C_new = solve_tridiagonal(K_D_diag, K_L_subdiag, rhs)

    # Stop if the prescribed flux produces a negative nodal concentration.
    if np.min(C_new) < 0:
        raise ValueError(
            f"Negative concentration at t = {time:g} s: "
            "the prescribed flux has taken the solution outside its physical range."
        )

    dC_dt = (C_new - C) / dt
    C = C_new

    surface_concentrations[step + 1] = C[-1, 0]
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
    f"R = {R:g} m | D = {D:g} m^2/s | c_initial = {c_initial:g} mol/m^3\n"
    f"F = {F:.8g} C/mol | elements = {N1} | steps = {n_steps}\n"
    f"dt = {dt:g} s | end time = {n_steps * dt:g} s\n"
    + textwrap.fill(f"i(t) = {current_expression} [A/m^2] (positive = insertion)", width=85)
    + f"\nj_const = {j_const:g} mol/(m^2*s) | i = {F * j_const:.8g} A/m^2"
)


def add_input_summary(figure):
    """Reserve a footer below the axes for the active simulation inputs."""
    footer_height = 0.25 + 0.18 * len(input_summary.splitlines())
    figure.set_size_inches(figure.get_figwidth(), figure.get_figheight() + footer_height)
    figure.text(
        0.02, 0.025, input_summary, ha="left", va="bottom", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="whitesmoke", edgecolor="0.8"),
    )
    figure.tight_layout(rect=(0, footer_height / figure.get_figheight(), 1, 1))


# Compare the profiles above and their signed difference below.
# The lower panel makes errors visible even when the profiles overlap.
fig, (ax, residual_ax) = plt.subplots(
    2, 1, sharex=True, figsize=(8, 8),
    gridspec_kw={"height_ratios": [1.2, 1]},
)
largest_error = 0.0
largest_normalised_error = None

for profile_index, (time, concentration) in enumerate(profiles):

    # Your numerical FEM solution.
    fem_values = concentration[:, 0]

    # Crank's solution at exactly the same radius nodes and time.
    analytical_values = crank_constant_flux(
        r=r,
        time=time,
        R=R,
        D=D,
        c0=c_initial,
        j=j_const,
        beta_roots=beta_roots
    )

    # Draw Crank as a line and FEM as sparse open circles of the same colour.
    colour = plt.get_cmap("tab10")(profile_index % 10)
    line, = ax.plot(r, analytical_values, linewidth=1.5, color=colour)
    ax.plot(
        r, fem_values, linestyle="none", marker="o", markersize=4,
        markerfacecolor="none", markeredgecolor=line.get_color(), color=colour,
        markevery=max(1, len(r) // 20),
    )

    residual = fem_values - analytical_values
    # Use one concentration-change scale for the WHOLE profile at each time.
    # Do not divide by local concentration, which can vanish near the centre.
    scale = np.max(np.abs(analytical_values - c_initial))
    if scale > 0:
        normalised_error = 100.0 * residual / scale
        max_normalised_error = np.max(np.abs(normalised_error))
        largest_normalised_error = max(
            largest_normalised_error or 0.0, max_normalised_error
        )
        normalised_label = f"{max_normalised_error:.6g}%"
        time_label = f"t = {time:g}"
    else:
        # Relative error is undefined at t=0 (or for zero applied current).
        normalised_error = np.full_like(residual, np.nan)
        normalised_label = "undefined (zero analytical change)"
        time_label = f"t = {time:g} (undefined)"
    residual_ax.plot(r, normalised_error, color=line.get_color(), label=time_label)
    max_error = np.max(np.abs(residual))
    largest_error = max(largest_error, max_error)

    print(
        f"t = {time:g}: "
        f"maximum absolute error = {max_error:.6e} mol/m^3; "
        f"maximum normalised error = {normalised_label}"
    )

# Method legend above; time legend below applies to both panels.
ax.plot([], [], "k-", label="Crank analytical solution")
ax.plot([], [], linestyle="none", marker="o", color="black",
        markerfacecolor="none", label="FEM (sampled nodes)")

ax.set_ylabel("Concentration [mol/m^3]")
ax.set_title("Constant surface flux: concentration and comparison error")
ax.ticklabel_format(axis="y", style="plain", useOffset=False)
ax.grid(True, alpha=0.3)
ax.legend()

residual_ax.axhline(0.0, color="black", linewidth=0.8, linestyle="--")
residual_ax.set_xlabel("Radius r [m]")
residual_ax.set_ylabel("Normalised error [%]")
normalised_summary = (
    f"Largest magnitude: {largest_normalised_error:.3g}%"
    if largest_normalised_error is not None
    else "Undefined: analytical concentration change is zero"
)
residual_ax.set_title(
    "100 * (FEM - Crank) / max_r |Crank - c_initial|\n"
    + normalised_summary + " (t = 0 undefined)",
    fontsize=10,
)
residual_ax.ticklabel_format(axis="y", style="plain", useOffset=False)
residual_ax.grid(True, alpha=0.3)
residual_ax.legend(title="Time (colours in both panels)", ncol=3, fontsize=8)
add_input_summary(fig)

# Save beside this script; also used to name the volume-average plot below.
plot_path = Path(__file__).with_name(plot_filename)
fig.savefig(plot_path, dpi=plot_dpi)
print(f"Saved concentration comparison to {plot_path}")

# Plot the volume average at every time step, including a current-balance reference.
average_fig, average_ax = plt.subplots(figsize=figure_size)
average_ax.plot(times, average_concentrations, label="FEM volume average")
average_ax.plot(
    times, current_balance_averages, "--", label="Current balance (backward Euler)"
)
# Analytical conservation law for uniform initial concentration and constant flux.
analytical_averages = c_initial + 3.0 * j_const * times / R
average_ax.plot(times, analytical_averages, ":", color="black",
                label="Analytical average: c0 + 3*j*t/R")
average_ax.set_xlabel("Time t [s]")
average_ax.set_ylabel("Volume-average concentration [mol/m^3]")
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

# Surface concentration c(R,t): at x = 1, every Crank spatial mode equals 1.
# The spatial polynomial becomes 1/2 - 3/10 = 1/5.
# Early times need more series terms. Choose enough that omitted modes have
# decay factors below about exp(-36) at the first positive time, dt.
surface_number_of_terms = max(
    number_of_terms, int(np.ceil(6 * R / (np.pi * np.sqrt(D * dt))))
)
surface_roots = (
    beta_roots if surface_number_of_terms == number_of_terms
    else calculate_crank_roots(surface_number_of_terms)
)
surface_series = np.zeros_like(times)
for beta in surface_roots:
    surface_series += np.exp(-D * beta**2 * times / R**2) / beta**2

analytical_surface = (
    c_initial + 3 * j_const * times / R
    + (j_const * R / D) * (1 / 5 - 2 * surface_series)
)
# Use the exact initial value rather than a finite-series approximation at t=0.
analytical_surface[0] = c_initial
surface_error = surface_concentrations - analytical_surface

surface_fig, (surface_ax, surface_error_ax) = plt.subplots(
    2, 1, sharex=True, figsize=(8, 7),
    gridspec_kw={"height_ratios": [1.2, 1]},
)
surface_ax.plot(times, analytical_surface, color="black",
                label="Crank analytical solution")
surface_ax.plot(
    times, surface_concentrations, linestyle="none", marker="o",
    markersize=4, markerfacecolor="none", color="tab:blue",
    markevery=max(1, n_steps // 25), label="FEM (sampled time steps)",
)
surface_ax.set_ylabel("Surface concentration [mol/m^3]")
surface_ax.set_title("Particle surface: concentration over time")
surface_ax.ticklabel_format(axis="y", style="plain", useOffset=False)
surface_ax.grid(True, alpha=0.3)
surface_ax.legend()

surface_error_ax.plot(times, surface_error, color="tab:blue")
surface_error_ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
surface_error_ax.set_xlabel("Time t [s]")
surface_error_ax.set_ylabel("Surface FEM - Crank [mol/m^3]")
surface_error_ax.set_title("Negative values mean FEM predicts a lower surface concentration",
                          fontsize=10)
surface_error_ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0),
                                 useOffset=False)
surface_error_ax.grid(True, alpha=0.3)
add_input_summary(surface_fig)

surface_plot_path = plot_path.with_name(f"{plot_path.stem}_surface.png")
surface_fig.savefig(surface_plot_path, dpi=plot_dpi)
surface_csv_path = plot_path.with_name(f"{plot_path.stem}_surface.csv")
np.savetxt(
    surface_csv_path,
    np.column_stack((times, surface_concentrations, analytical_surface, surface_error)),
    delimiter=",", comments="",
    header="time_s,FEM_surface_mol_per_m3,Crank_surface_mol_per_m3,FEM_minus_Crank_mol_per_m3",
)
print(f"Saved surface plot to {surface_plot_path}")
print(f"Saved surface data to {surface_csv_path}")
print(f"Surface reference terms: {surface_number_of_terms}")
print(f"Final FEM surface concentration: {surface_concentrations[-1]:.10g}")
print(f"Final Crank surface concentration: {analytical_surface[-1]:.10g}")
print(f"Maximum absolute surface error over all time steps: {np.max(np.abs(surface_error)):.3e}")
plt.show()
