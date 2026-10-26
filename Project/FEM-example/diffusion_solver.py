"""Linear spherical FEM with constant D and current density."""
import numpy as np

FARADAY = 96485.3329


def solve_diffusion(n_elements, radius, diffusivity, current, initial, dt, end_time):
    """Return nodes, times and concentration changes using backward Euler.

    The consistent mass and diffusion matrices are stored as three diagonals.
    Solving for change from the uniform initial value avoids cancellation when
    the imposed current causes very small concentration changes.
    """
    n_elements = int(n_elements)
    r = np.linspace(0.0, radius, n_elements + 1)
    h = radius / n_elements
    left, right = r[:-1], r[1:]
    factor = 4 * np.pi * h / 60
    m11 = factor * (12 * left**2 + 6 * left * right + 2 * right**2)
    m12 = factor * (3 * left**2 + 4 * left * right + 3 * right**2)
    m22 = factor * (2 * left**2 + 6 * left * right + 12 * right**2)
    a = 4 * np.pi * diffusivity * (right**3 - left**3) / (3 * h**2)

    mass_diag = np.zeros(len(r))
    mass_diag[:-1] += m11
    mass_diag[1:] += m22
    stiff_diag = np.zeros(len(r))
    stiff_diag[:-1] += a
    stiff_diag[1:] += a

    # ceil ensures the requested time step is never exceeded. Shorten the last
    # step, if necessary, to land exactly on the requested end time.
    count = max(1, int(np.ceil(end_time / dt)))
    times = np.minimum(np.arange(count + 1) * dt, end_time)
    times[-1] = end_time
    changes = np.zeros((count + 1, len(r)))
    boundary_load = 4 * np.pi * radius**2 * current / FARADAY
    factors = {}

    for step in range(count):
        step_dt = float(times[step + 1] - times[step])
        # Reuse a tridiagonal factorisation when the time step is unchanged.
        key = round(step_dt / dt, 12)
        if key not in factors:
            off = m12 - step_dt * a
            pivot = mass_diag + step_dt * stiff_diag
            upper = np.empty(n_elements)
            for k in range(n_elements):
                upper[k] = off[k] / pivot[k]
                pivot[k + 1] -= off[k] * upper[k]
            factors[key] = (off, pivot, upper)
        off, pivot, upper = factors[key]

        old = changes[step]
        rhs = mass_diag * old
        rhs[:-1] += m12 * old[1:]
        rhs[1:] += m12 * old[:-1]
        rhs[-1] += step_dt * boundary_load
        new = np.empty(len(r))
        new[0] = rhs[0] / pivot[0]
        for k in range(1, len(r)):
            new[k] = (rhs[k] - off[k - 1] * new[k - 1]) / pivot[k]
        for k in range(n_elements - 1, -1, -1):
            new[k] -= upper[k] * new[k + 1]
        changes[step + 1] = new

    weights = mass_diag.copy()
    weights[:-1] += m12
    weights[1:] += m12
    average_change = changes @ weights / (4 * np.pi * radius**3 / 3)
    return r, times, changes, average_change

