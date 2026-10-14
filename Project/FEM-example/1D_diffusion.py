import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

N1 = 20          # Number of elements
N0 = N1 + 1       # Number of vertices
R=1.0
h = R / N1      # Element length
pi = np.pi
D = 1.0
r = np.linspace(0.0, R, N0)
F = 96485.3329
# Each element connects two neighbouring vertices
elements = np.column_stack((
    np.arange(N1),
    np.arange(1, N0),
))
c_initial = 1.0
C = np.full((len(r), 1), c_initial)

A = np.zeros((N0, N0))  # Global matrix: N0 × N0
B = np.zeros((N0, 1))  # Global vector: N0 × 1
M = np.zeros((N0, N0))  # Global mass matrix: N0 × N0
def local_mass_matrix(rL, rR):
    h = rR - rL

    M_e = (4 * np.pi * h / 60) * np.array([
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
        4 * np.pi * D * (r_right**3 - r_left**3) / (3 * h**2)
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


# Prescribed current density. Positive current adds lithium to the particle.
def i(t):
    return 1.0  # Set your current density here, or use an expression involving t.

# Time stepping
dt = 0.01
n_steps = 100
K = M + dt * A

# Store the initial profile and four later profiles for plotting.
snapshot_steps = set(np.linspace(0, n_steps, 5, dtype=int))
profiles = [(0.0, C.copy())]

for step in range(n_steps):
    time = (step + 1) * dt
    # Backward Euler evaluates the surface current at the new time.
    B[-1, 0] = 4 * np.pi * R**2 * i(time) / F
    rhs = M @ C + dt * B
    C_new = np.linalg.solve(K, rhs)

    dC_dt = (C_new - C) / dt
    C = C_new

    if step + 1 in snapshot_steps:
        profiles.append(((step + 1) * dt, C.copy()))



# Plot concentration against radius at the saved times.
fig, ax = plt.subplots(figsize=(7, 4.5))
for time, concentration in profiles:
    ax.plot(r, concentration[:, 0], label=f"t = {time:.2f}")

ax.set_xlabel("Radius r")
ax.set_ylabel("Concentration C")
ax.set_title("Spherical diffusion: concentration profiles")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()

plot_path = Path(__file__).with_name("concentration_profiles.png")
fig.savefig(plot_path, dpi=180)
print(f"Saved plot to {plot_path}")
plt.show()
