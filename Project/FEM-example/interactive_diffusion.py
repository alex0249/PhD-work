"""Interactive linear FEM for spherical diffusion.

Run with the same Python environment as 1D_diffusion.py.
Uses NumPy and Matplotlib only. Positive surface current density inserts lithium.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, CheckButtons, RadioButtons, Slider, TextBox
from matplotlib.ticker import ScalarFormatter

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


def plain_ticks(ax):
    """Show full concentration values, without hidden additive offsets."""
    for axis in (ax.xaxis, ax.yaxis):
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        axis.set_major_formatter(formatter)


class DiffusionExplorer:
    def __init__(self):
        self.fig = plt.figure(figsize=(12, 9))
        self.fig.canvas.manager.set_window_title("Spherical diffusion explorer")
        self.fig.suptitle("Spherical diffusion | positive current enters the particle", y=0.99)
        self.fields = {}
        controls = [
            ("elements", "Elements", "1000"),
            ("radius", "Radius [m]", "1"),
            ("diffusivity", "D [m²/s]", "1"),
            ("current", "Current density [A/m²]", "1"),
            ("initial", "Initial c [mol/m³]", "1"),
            ("dt", "Time step [s]", "0.0001"),
            ("end", "End time [s]", "0.01"),
        ]
        for index, (name, label, default) in enumerate(controls):
            column, row = index % 4, index // 4
            xpos, ypos = 0.07 + column * 0.235, 0.91 - row * 0.085
            self.fig.text(xpos, ypos + 0.034, label, fontsize=9)
            textbox = TextBox(self.fig.add_axes([xpos, ypos, 0.19, 0.029]), "", initial=default)
            self.fields[name] = textbox

        self.run_button = Button(self.fig.add_axes([0.775, 0.825, 0.19, 0.033]), "Run simulation")
        self.run_button.on_clicked(self.run)
        self.status = self.fig.text(0.07, 0.79, "", fontsize=9)
        self.profile_ax = self.fig.add_axes([0.1, 0.455, 0.86, 0.285])
        self.history_ax = self.fig.add_axes([0.1, 0.19, 0.86, 0.175])
        self.slider_ax = self.fig.add_axes([0.16, 0.108, 0.70, 0.025])
        self.slider = None
        self.slider_cid = None
        self.mode = RadioButtons(self.fig.add_axes([0.07, 0.01, 0.28, 0.072]),
                                 ["Change from initial", "Absolute concentration"], active=0)
        self.options = CheckButtons(self.fig.add_axes([0.38, 0.01, 0.29, 0.072]),
                                    ["Overlay early profiles", "Zoom to outer 20%"], [True, False])
        self.mode.on_clicked(self.redraw)
        self.options.on_clicked(self.redraw)
        self.time_label = self.fig.text(0.73, 0.047, "", fontsize=10)
        self.result = None
        self.run()

    def run(self, _event=None):
        try:
            values = {key: float(box.text) for key, box in self.fields.items()}
            if not all(np.isfinite(v) for v in values.values()):
                raise ValueError("Inputs must be finite numbers.")
            n = int(values["elements"])
            if n != values["elements"] or not 2 <= n <= 5000:
                raise ValueError("Elements must be an integer between 2 and 5000.")
            if min(values["radius"], values["diffusivity"], values["dt"], values["end"]) <= 0:
                raise ValueError("Radius, D, time step and end time must be positive.")
            if values["initial"] < 0:
                raise ValueError("Initial concentration must be nonnegative.")
            steps = max(1, int(np.ceil(values["end"] / values["dt"])))
            if steps > 20000 or (steps + 1) * (n + 1) > 3_000_000:
                raise ValueError("Use fewer elements/steps (limit: 20,000 steps and 3 million stored values).")
            self.status.set_text("Calculating…")
            self.status.set_color("black")
            self.fig.canvas.draw_idle()
            self.fig.canvas.flush_events()
            result = solve_diffusion(n, values["radius"], values["diffusivity"],
                                     values["current"], values["initial"], values["dt"], values["end"])
            self.result, self.values = result, values
            if self.slider is not None:
                self.slider.disconnect(self.slider_cid)
                self.slider.disconnect_events()
            self.slider_ax.clear()
            self.slider = Slider(self.slider_ax, "Time step", 0, steps, valinit=steps, valstep=1)
            self.slider_cid = self.slider.on_changed(self.redraw)
            self.redraw()
        except (ValueError, OverflowError) as error:
            self.status.set_text(str(error))
            self.status.set_color("firebrick")
            self.fig.canvas.draw_idle()

    def redraw(self, _event=None):
        if self.result is None:
            return
        r, times, changes, averages = self.result
        selected = int(self.slider.val)
        absolute = self.mode.value_selected == "Absolute concentration"
        overlay, zoom = self.options.get_status()
        initial = self.values["initial"]
        offset = initial if absolute else 0.0
        ylabel = "Concentration [mol/m³]" if absolute else "Concentration change [mol/m³]"
        self.profile_ax.clear()
        self.history_ax.clear()
        if overlay and selected > 0:
            candidates = [0] + np.geomspace(1, selected, min(5, selected)).round().astype(int).tolist()
            for index in sorted(set(candidates)):
                if index != selected:
                    self.profile_ax.plot(r / r[-1], changes[index] + offset,
                                         linewidth=1.2, alpha=0.7,
                                         label=f"t = {times[index]:.6g} s")
        self.profile_ax.plot(r / r[-1], changes[selected] + offset,
                             linewidth=2.4, label=f"Selected: {times[selected]:.6g} s")
        self.profile_ax.set_xlim(0.8 if zoom else 0.0, 1.0)
        self.profile_ax.set_xlabel("Fraction of radius, r/R")
        self.profile_ax.set_ylabel(ylabel)
        self.profile_ax.set_title("Radial profile", loc="left", fontsize=11)
        self.profile_ax.legend(fontsize=8, loc="best", ncol=2)

        self.history_ax.plot(times, changes[:, -1] + offset, label="Surface", linewidth=1.8)
        self.history_ax.plot(times, changes[:, 0] + offset, label="Centre", linewidth=1.4)
        self.history_ax.plot(times, averages + offset, "--", label="Volume average", linewidth=1.4)
        self.history_ax.axvline(times[selected], color="0.5", linewidth=1)
        self.history_ax.set_xlabel("Time [s]")
        self.history_ax.set_ylabel(ylabel)
        self.history_ax.set_xlim(0, times[-1])
        self.history_ax.legend(fontsize=8, loc="best", ncol=3)
        for ax in (self.profile_ax, self.history_ax):
            ax.grid(True, alpha=0.2)
            plain_ticks(ax)
        expected = 3 * self.values["current"] * times[-1] / (FARADAY * r[-1])
        balance_error = abs(averages[-1] - expected)
        self.status.set_text(
            f"{len(r)-1} elements · {len(times)-1} steps · "
            f"Final average change {averages[-1]:.6g} · Balance error {balance_error:.2g}"
        )
        self.status.set_color("black")
        self.time_label.set_text(f"t = {times[selected]:.6g} s\nSurface c = {initial + changes[selected, -1]:.10g}")
        self.fig.canvas.draw_idle()


if __name__ == "__main__":
    explorer = DiffusionExplorer()
    plt.show()
