import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ==================== DIMENSIONLESS INPUTS ====================
# c_star = c/c_max (zero reference concentration).
# I_star = R*i/(F*D*c_max) = R*j/(D*c_max).
# Positive I_star inserts lithium. The Crank comparison assumes constant current.
I_star = 10.0              # Main physical control parameter, dimensionless.
c_initial_star = 0.0      # Uniform initial concentration c_initial/c_max.
# Scaling does not enforce c_star <= 1; this linear model has no saturation law.

# Numerical settings (also dimensionless).
N1 = 1000                # Number of radial elements.
dt_star = 1e-5           # dt_star = D*dt/R**2.
n_steps = 10000          # End time t_star = n_steps*dt_star.
snapshot_steps = {0, 100, 500, 1000, 2000, 5000, n_steps}
plot_dpi = 180
plot_filename = "non_dimensional_comparison.png"
# ==================== END OF INPUTS ====================

# Validate inputs
if not np.isfinite(I_star):
    raise ValueError("I_star must be finite.")
if not np.isfinite(c_initial_star) or c_initial_star < 0:
    raise ValueError("c_initial_star must be finite and nonnegative.")
if not np.isfinite(dt_star) or dt_star <= 0:
    raise ValueError("dt_star must be finite and positive.")
if not isinstance(N1, (int, np.integer)) or N1 < 1:
    raise ValueError("N1 must be a positive integer.")
if not isinstance(n_steps, (int, np.integer)) or n_steps < 1:
    raise ValueError("n_steps must be a positive integer.")

# Set non-dimensional current input function
def current_star(t_star):
    """Constant dimensionless inward current, shared by FEM and Crank."""
    return I_star
## Time dependent current
# def current_star(t_star):
#     return 3 * I_star * t_star

# Derived dimensionless mesh quantities
N0 = N1 + 1       # Number of vertices
h_star = 1.0 / N1  # Dimensionless element length
pi = np.pi
r_star = np.linspace(0.0, 1.0, N0)
# Each element connects two neighbouring vertices
elements = np.column_stack((
    np.arange(N1),
    np.arange(1, N0),
))
# Initialise the FEM solution vector with the uniform initial concentration.
C_star = np.full((N0, 1), c_initial_star)

# Symmetric tridiagonal matrices: store the diagonal and one off-diagonal.
A_diag_star = np.zeros(N0)
A_off_star = np.zeros(N1)
M_diag_star = np.zeros(N0)
M_off_star = np.zeros(N1)
B_star = np.zeros((N0, 1))
# Three-point Gauss-Legendre quadrature on the reference element [-1, 1].
gauss_points = np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)])
gauss_weights = np.array([5 / 9, 8 / 9, 5 / 9])


def local_element_matrices(rL_star, rR_star):
    """Integrate A_e and M_e numerically, including the radial r*^2 weight."""
    # Define the element length
    element_h_star = rR_star - rL_star
    # Define the Jacboian for the transformation from reference element [-1, 1] to [rL_star, rR_star]
    jacobian = element_h_star / 2  # dr*/dxi
    # Result for the derivative of the shape functions (phi = (1 +/- xi)/2) with respect to the reference coordinate xi
    dphi_dxi = np.array([-0.5, 0.5])
    # Transform the derivative to the physical coordinate r* using the chain rule: dphi/dr* = dphi/dxi * dxi/dr* = dphi/dxi / (dr*/dxi)
    dphi_dr_star = dphi_dxi / jacobian
    # Initialise the local matrices for the element
    A_e = np.zeros((2, 2))
    M_e = np.zeros((2, 2))
    # Perform numerical integration using Gauss-Legendre quadrature
    for xi_q, w_q in zip(gauss_points, gauss_weights):
        # Linear shape functions at this quadrature point.
        # Evaluate the shape functions at the parent-coordinate Gauss point xi_q
        phi = np.array([(1 - xi_q) / 2, (1 + xi_q) / 2])
        # Find corresponding physical radial coordinate
        r_q_star = rL_star + (1 + xi_q) * jacobian
        # Find integration weight for this quadrature point
        integration_weight = w_q * jacobian * r_q_star**2

        # M_ij += w_q (h*/2) phi_i phi_j (r_q*)^2.
        M_e += integration_weight * np.outer(phi, phi)
        # A_ij += w_q (h*/2) (dphi_i/dr*) (dphi_j/dr*) (r_q*)^2.
        A_e += integration_weight * np.outer(dphi_dr_star, dphi_dr_star)

    return A_e, M_e

# Asssemble the global matrices A and M from the local element matrices.
for e in range(N1):
    # Get the left and right radial coordinates of the element
    r_left = r_star[e]
    r_right = r_star[e + 1]

    # Compute the local element matrices for this element
    A_e, M_e = local_element_matrices(r_left, r_right)

    # Assemble directly into the nonzero diagonals.
    A_diag_star[e] += A_e[0, 0]
    A_diag_star[e + 1] += A_e[1, 1]
    A_off_star[e] += A_e[0, 1]
    M_diag_star[e] += M_e[0, 0]
    M_diag_star[e + 1] += M_e[1, 1]
    M_off_star[e] += M_e[0, 1]


def tridiagonal_matvec(diagonal, off_diagonal, values):
    """Multiply a symmetric tridiagonal matrix by a column vector."""
    # Multiply each row by the diagonal and off-diagonal contributions.
    result = diagonal[:, None] * values  
    result[:-1] += off_diagonal[:, None] * values[1:]
    result[1:] += off_diagonal[:, None] * values[:-1]
    return result


def factor_tridiagonal(diagonal, off_diagonal):
    """Compute K = L D L.T for symmetric positive-definite K.

    D_diag stores the diagonal of D; L_subdiag stores the subdiagonal of L.
    L has ones on its diagonal. Both factors are stored as 1D arrays.
    """
    # Copy the diagonal to avoid modifying the input.
    D_diag = diagonal.copy()
    # Initialise the subdiagonal of L.
    L_subdiag = np.empty_like(off_diagonal)
    # Loop over each row to compute the K = L D L.T factorisation.
    for row in range(len(off_diagonal)):
        # Compute the subdiagonal of L and update the diagonal of D.
        L_subdiag[row] = off_diagonal[row] / D_diag[row]
        D_diag[row + 1] -= L_subdiag[row] * off_diagonal[row]
    # Check that the resulting diagonal is positive definite.
    if not np.all(np.isfinite(D_diag)) or np.any(D_diag <= 0):
        raise ValueError("The diffusion matrix must be positive definite.")
    return D_diag, L_subdiag


def solve_tridiagonal(D_diag, L_subdiag, rhs):
    """Solve K x = rhs using D_diag and L_subdiag from K = L D L.T."""
    # Copy the right-hand side to avoid modifying the input.
    solution = rhs[:, 0].copy()
    # Forward substitution: L y = rhs. Lower triangular system here
    for row in range(1, len(D_diag)):
        # Starting at the second row, subtract the previous, already 
        # updated value multiplied by the corresponding entry of L
        solution[row] -= L_subdiag[row - 1] * solution[row - 1]
    # Diagonal solve: D z = y. 
    # Takes the y from the forward sub step and finds z.
    # Divide each value by the 
    # corresponding diagonal entry of D
    solution /= D_diag
    # Back substitution: L.T x = z. 
    # Takes the z from the diagonal solve and finds x.
    for row in range(len(D_diag) - 2, -1, -1):
        solution[row] -= L_subdiag[row] * solution[row + 1]
    # Solution to Ax = rhs ( or MC_dot + AC = B specifically) is now in solution. Return as a column vector.
    return solution[:, None]

def calculate_crank_roots(number_of_terms):
    """Find the positive roots beta_n of tan(beta_n) = beta_n.
    Uses bisection to find the roots in the intervals (n*pi, n*pi + pi/2) for n = 1, 2, ..., number_of_terms."""

    roots = np.zeros(number_of_terms)
    # Loop over each term to find the root.
    for n in range(1, number_of_terms + 1):

        # The nth positive root lies in this interval.
        left = n * np.pi
        right = n * np.pi + np.pi / 2

        # Repeatedly halve the interval containing the root.
        for iteration in range(60):
            # Compute the midpoint of the current interval.
            middle = (left + right) / 2

            # Evaluate the function at the left and middle points 
            # to determine which side of the interval contains the root.
            # Function is tan(beta) = beta so
            # Evaluate the function at the left point
            value_left = np.sin(left) - left * np.cos(left)
            # Evaluate the function at the middle point
            value_middle = np.sin(middle) - middle * np.cos(middle)
            # If the function has the same sign at left and middle, 
            # the root is in the right half of the interval.
            if value_left * value_middle > 0:
                left = middle
            else:
                # Otherwise, the root is in the left half of the interval.
                right = middle
        # The root is now isolated in the interval [left, right].
        roots[n - 1] = (left + right) / 2
    return roots
number_of_terms = 80
# Very early times require more modes to cancel the spatial polynomial accurately.
# Use enough terms that omitted modes decay below about exp(-36) at dt_star.
number_of_terms = max(
    number_of_terms, int(np.ceil(6 / (np.pi * np.sqrt(dt_star))))
) # The six here is just selected 
beta_roots = calculate_crank_roots(number_of_terms)
print(f"Crank reference terms: {number_of_terms}")

def crank_constant_flux_star(r_star, t_star, c0_star, I_star, beta_roots):
    """Analytical concentration for a sphere with constant inward flux.

    c*(r*, t*) = c0* + I* [3 t* + r*^2/2 - 3/10
                    - 2 sum_n g_n(r*) exp(-beta_n^2 t*) / beta_n^2]
    where g_n(r*) = sin(beta_n r*) / (r* sin(beta_n)) and
    tan(beta_n) = beta_n. The centre uses the finite limit of g_n.

    This solves dc*/dt* = (1/r*^2) d/dr*(r*^2 dc*/dr*) with
    dc*/dr*(0, t*) = 0, dc*/dr*(1, t*) = I*, and c*(r*, 0) = c0*.
    The infinite series is approximated using the supplied beta_roots.
    """

    x = np.asarray(r_star, dtype=float)

    # Check whether evaluating at the initial time.
    if t_star == 0:
        # Return the uniform initial concentration for all radii.
        return np.full_like(x, c0_star)

    # Accumulate the transient series, one mode at a time.
    series = np.zeros_like(x)

    for beta in beta_roots:
        # Each spatial mode is gn_x used to build the concentration profile. The series is summed over all modes.
        # g_n(x) = sin(beta*x) / (x*sin(beta)) = sinc(beta*x/pi) / sinc(beta/pi)
        # Using sinc also handles the centre x = 0 safely as np.sinc(0) = 1.
        spatial_mode = (
            np.sinc(beta * x / np.pi)
            / np.sinc(beta / np.pi)
        )
        # Inverse laplacian of the non-zero poles
        decay = np.exp(-beta**2 * t_star)
        # Add the contribution of this mode to the series.
        series += spatial_mode * decay / beta**2

    # The volume-average concentration increases with insertion.
    mean_concentration = c0_star + 3 * I_star * t_star

    # Spatial variation relative to the volume average.
    spatial_correction = I_star * (
        x**2 / 2 - 3 / 10 - 2 * series
    )

    return mean_concentration + spatial_correction

# K_star is constant: factor it once and reuse at every step.
K_diag_star = M_diag_star + dt_star * A_diag_star
K_off_star = M_off_star + dt_star * A_off_star
K_D_diag_star, K_L_subdiag_star = factor_tridiagonal(K_diag_star, K_off_star)

# Store the initial profile; later profiles_star are saved at snapshot_steps.
profiles_star = [(0.0, C_star.copy())]

# Array for time series
times_star = np.arange(n_steps + 1) * dt_star
# Record the outermost FEM node at EVERY step, not just the profile snapshots.
surface_concentrations_star = np.empty(n_steps + 1)
surface_concentrations_star[0] = C_star[-1, 0]

# Time stepping loop
for step in range(n_steps):
    # Advance time by one step
    time_star = (step + 1) * dt_star
    # Backward Euler evaluates the applied current at end of the last time step, at the new time_star.
    B_star[-1, 0] = current_star(time_star) # Last row of B_star is the applied flux at the surface
    rhs = tridiagonal_matvec(M_diag_star, M_off_star, C_star) + dt_star * B_star
    C_new_star = solve_tridiagonal(K_D_diag_star, K_L_subdiag_star, rhs)

    # Stop if the prescribed flux produces a negative nodal concentration.
    if np.min(C_new_star) < 0:
        raise ValueError(
            f"Negative dimensionless concentration at t_star = {time_star:g}: "
            "the prescribed flux has taken the solution outside its physical range."
        )
    # Calculate the rate of change of the concentration.
    dC_star_dt_star = (C_new_star - C_star) / dt_star
    # Update the concentration for the next time step.
    C_star = C_new_star
    # Record the surface concentration at the outermost node.
    surface_concentrations_star[step + 1] = C_star[-1, 0]
    # Save the profile at the requested snapshot steps or at the final step.
    if step + 1 in snapshot_steps or step + 1 == n_steps:
        profiles_star.append(((step + 1) * dt_star, C_star.copy()))



# ==================== PLOTTING ====================

# All plots and exported data below are dimensionless.
input_summary = (
    f"I_star = {I_star:g} | c_initial_star = {c_initial_star:g}\n"
    f"elements = {N1} | steps = {n_steps} | dt_star = {dt_star:g}\n"
    f"end t_star = {n_steps * dt_star:g} | positive current = insertion"
)

def add_input_summary(figure):
    footer_height = 0.25 + 0.18 * len(input_summary.splitlines())
    figure.set_size_inches(figure.get_figwidth(), figure.get_figheight() + footer_height)
    figure.text(0.02, 0.025, input_summary, ha="left", va="bottom", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.4", facecolor="whitesmoke", edgecolor="0.8"))
    figure.tight_layout(rect=(0, footer_height / figure.get_figheight(), 1, 1))

fig, (ax, residual_ax) = plt.subplots(
    2, 1, sharex=True, figsize=(8, 8), gridspec_kw={"height_ratios": [1.2, 1]})
largest_normalised_error = None
profile_comparison_rows = []
for profile_index, (t_star, concentration_star) in enumerate(profiles_star):
    fem_values = concentration_star[:, 0]
    analytical_values = crank_constant_flux_star(
        r_star, t_star, c_initial_star, I_star, beta_roots)
    colour = plt.get_cmap("tab10")(profile_index % 10)
    ax.plot(r_star, analytical_values, linewidth=1.5, color=colour)
    ax.plot(r_star, fem_values, linestyle="none", marker="o", markersize=4,
            markerfacecolor="none", color=colour, markevery=max(1, N0 // 20))
    residual = fem_values - analytical_values
    # Keep values at every mesh node for direct numerical comparison.
    profile_comparison_rows.append(np.column_stack((
        np.full_like(r_star, t_star), r_star,
        fem_values, analytical_values, residual,
    )))
    scale = np.max(np.abs(analytical_values - c_initial_star))
    if scale > 0:
        normalised_error = 100 * residual / scale
        error_percent = np.max(np.abs(normalised_error))
        largest_normalised_error = max(largest_normalised_error or 0., error_percent)
        error_label = f"{error_percent:.6g}%"
    else:
        normalised_error = np.full_like(residual, np.nan)
        error_label = "undefined (zero analytical change)"
    residual_ax.plot(r_star, normalised_error, color=colour, label=f"t* = {t_star:g}")
    print(f"t_star = {t_star:g}: max |FEM - Crank| = {np.max(np.abs(residual)):.6e} "
          f"(dimensionless); maximum normalised error = {error_label}")

ax.plot([], [], "k-", label="Crank analytical solution")
ax.plot([], [], linestyle="none", marker="o", color="black", markerfacecolor="none",
        label="FEM (sampled nodes)")
ax.set_ylabel("Concentration c* = c/c_max")
ax.set_title("Dimensionless constant-flux diffusion")
ax.ticklabel_format(axis="y", style="plain", useOffset=False)
ax.grid(True, alpha=0.3)
ax.legend()
residual_ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
residual_ax.set_xlabel("Radius r* = r/R")
residual_ax.set_ylabel("Normalised error [%]")
error_summary = (f"Largest magnitude: {largest_normalised_error:.3g}%"
                 if largest_normalised_error is not None else "Undefined: zero analytical change")
residual_ax.set_title("100 * (FEM - Crank) / max_r* |Crank - c_initial*|\n" + error_summary,
                      fontsize=10)
residual_ax.grid(True, alpha=0.3)
residual_ax.legend(title="Dimensionless time", ncol=3, fontsize=8)
add_input_summary(fig)
plot_path = Path(__file__).with_name(plot_filename)
fig.savefig(plot_path, dpi=plot_dpi)
profile_csv_path = plot_path.with_name(f"{plot_path.stem}_profiles.csv")
np.savetxt(
    profile_csv_path, np.vstack(profile_comparison_rows),
    delimiter=",", comments="",
    header="time_star,radius_star,FEM_c_star,Crank_c_star,FEM_minus_Crank_c_star",
)
# At r_star=1 every Crank spatial mode is one.
surface_series = np.zeros_like(times_star)
for beta in beta_roots:
    surface_series += np.exp(-beta**2 * times_star) / beta**2
analytical_surface_star = (c_initial_star + 3 * I_star * times_star
                           + I_star * (1/5 - 2*surface_series))
analytical_surface_star[0] = c_initial_star
surface_error_star = surface_concentrations_star - analytical_surface_star
surface_fig, (surface_ax, surface_error_ax) = plt.subplots(
    2, 1, sharex=True, figsize=(8, 7), gridspec_kw={"height_ratios": [1.2, 1]})
surface_ax.plot(times_star, analytical_surface_star, color="black", label="Crank analytical solution")
surface_ax.plot(times_star, surface_concentrations_star, linestyle="none", marker="o",
                markersize=4, markerfacecolor="none", color="tab:blue",
                markevery=max(1, n_steps // 25), label="FEM (sampled time steps)")
surface_ax.set_ylabel("Surface concentration c*(1,t*)")
surface_ax.set_title("Dimensionless particle surface concentration")
surface_ax.grid(True, alpha=0.3)
surface_ax.legend()
surface_error_ax.plot(times_star, surface_error_star, color="tab:blue")
surface_error_ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
surface_error_ax.set_xlabel("Time t* = Dt/R^2")
surface_error_ax.set_ylabel("Surface FEM - Crank (dimensionless)")
surface_error_ax.set_title("Negative values mean FEM predicts a lower surface concentration", fontsize=10)
surface_error_ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0), useOffset=False)
surface_error_ax.grid(True, alpha=0.3)
add_input_summary(surface_fig)
surface_plot_path = plot_path.with_name(f"{plot_path.stem}_surface.png")
surface_fig.savefig(surface_plot_path, dpi=plot_dpi)
surface_csv_path = plot_path.with_name(f"{plot_path.stem}_surface.csv")
np.savetxt(surface_csv_path,
           np.column_stack((times_star, surface_concentrations_star, analytical_surface_star, surface_error_star)),
           delimiter=",", comments="",
           header="time_star,FEM_surface_c_star,Crank_surface_c_star,FEM_minus_Crank_c_star")
for output_path in (plot_path, profile_csv_path,
                    surface_plot_path, surface_csv_path):
    print(f"Saved {output_path}")
print(f"Final FEM c_star at surface: {surface_concentrations_star[-1]:.10g}")
print(f"Final Crank c_star at surface: {analytical_surface_star[-1]:.10g}")
print(f"Maximum absolute dimensionless surface error: {np.max(np.abs(surface_error_star)):.3e}")
plt.show()
