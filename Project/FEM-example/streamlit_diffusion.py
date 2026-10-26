"""Run with: python -m streamlit run streamlit_diffusion.py"""
import io
import time
import numpy as np
import plotly.graph_objects as go
import streamlit as st
from diffusion_solver import solve_diffusion, FARADAY


st.set_page_config(page_title="Spherical diffusion", page_icon="◉", layout="wide")
st.title("Spherical diffusion")
st.caption("Edit the settings, then press Run simulation. Results remain unchanged until you run. Positive current density adds lithium.")


@st.cache_data(max_entries=4, show_spinner=False)
def calculate(n, radius, diffusivity, current, initial, dt, end):
    return solve_diffusion(n, radius, diffusivity, current, initial, dt, end)


previous = st.session_state.get("run_settings", {})
horizon = float(previous.get("end", 0.01))
last_view = st.session_state.get("view_settings", {})
with st.form("simulation_controls", enter_to_submit=False):
    columns = st.columns(4)
    with columns[0]:
        n = st.number_input("Elements", min_value=2, max_value=5000, value=1000, step=10, key="elements")
        initial = st.number_input("Initial concentration [mol/m³]", min_value=0.0, value=1.0, format="%.10g", key="initial")
    with columns[1]:
        radius = st.number_input("Radius [m]", min_value=0.0, value=1.0, format="%.10g", key="radius")
        dt = st.number_input("Time step Δt [s]", min_value=0.0, value=0.0001, format="%.10g", key="dt")
    with columns[2]:
        diffusivity = st.number_input("Diffusivity D [m²/s]", min_value=0.0, value=1.0, format="%.10g", key="D")
        end = st.number_input("End time [s]", min_value=0.0, value=0.01, format="%.10g", key="end")
    with columns[3]:
        current = st.number_input("Surface current density [A/m²]", value=1.0, format="%.10g", key="current")
        mode = st.selectbox("Vertical axis", ["Change from initial", "Absolute concentration"], key="mode")
    left, right = st.columns([2, 1])
    with left:
        requested_time = st.slider(
            "Display time [s]", min_value=0.0, max_value=horizon,
            value=min(float(last_view.get("time", horizon)), horizon),
            step=horizon / 1000, format="%.8f", key=f"display_time_{horizon:.12g}",
            help="Time in seconds, not a step number. The range follows the most recent run's end time."
        )
    with right:
        final_time = st.checkbox("Show final time", value=True, key="final_time")
        overlay = st.checkbox("Overlay early profiles", value=True, key="overlay")
        zoom = st.checkbox("Zoom to outer 20% of radius", value=False, key="zoom")
    submitted = st.form_submit_button("Run simulation", type="primary")

if submitted:
    settings = dict(n=int(n), radius=float(radius), diffusivity=float(diffusivity),
                    current=float(current), initial=float(initial), dt=float(dt), end=float(end))
    try:
        if not all(np.isfinite(v) for v in settings.values()):
            raise ValueError("Enter finite numbers for all inputs.")
        if min(radius, diffusivity, dt, end) <= 0:
            raise ValueError("Radius, diffusivity, time step and end time must be greater than zero.")
        steps = int(np.ceil(end / dt))
        if steps > 20000 or (steps + 1) * (n + 1) > 3_000_000:
            raise ValueError("Reduce the number of elements or steps: maximum 20,000 steps and 3 million stored values.")
        start = time.perf_counter()
        with st.spinner("Calculating concentration profiles…"):
            result = calculate(**settings)
        st.session_state["result"] = result
        st.session_state["run_settings"] = settings
        st.session_state["view_settings"] = dict(
            time=end if final_time else min(requested_time, end),
            mode=mode, overlay=overlay, zoom=zoom
        )
        st.session_state["elapsed"] = time.perf_counter() - start
        st.rerun()
    except (ValueError, OverflowError) as error:
        st.error(str(error))

if "result" not in st.session_state:
    st.info("Press Run simulation to display your first result.")
    st.stop()

# All charts use the last submitted settings, never the current draft widgets.
r, times, changes, averages = st.session_state["result"]
p = st.session_state["run_settings"]
view = st.session_state["view_settings"]
index = int(np.argmin(abs(times - view["time"])))
absolute = view["mode"] == "Absolute concentration"
offset = p["initial"] if absolute else 0.0
axis_title = "Concentration [mol/m³]" if absolute else "Concentration change [mol/m³]"
expected = 3 * p["current"] * times[-1] / (FARADAY * p["radius"])
metrics = st.columns(4)
metrics[0].metric("Displayed time [s]", f"{times[index]:.8g}")
metrics[1].metric("Surface concentration [mol/m³]", f"{p['initial'] + changes[index, -1]:.10g}")
metrics[2].metric("Volume-average change [mol/m³]", f"{averages[index]:.6g}")
metrics[3].metric("Final balance error [mol/m³]", f"{abs(averages[-1] - expected):.2g}")
st.caption(
    f"Showing the last submitted run: {p['n']} elements · {len(times)-1} steps · "
    f"Δt = {p['dt']:.8g} s · End time = {p['end']:.8g} s · "
    f"Calculation/cache retrieval: {st.session_state['elapsed']:.3f} s"
)

profile = go.Figure()
if view["overlay"] and index > 0:
    candidates = [0] + np.geomspace(1, index, min(5, index)).round().astype(int).tolist()
    for k in sorted(set(candidates)):
        if k != index:
            profile.add_trace(go.Scatter(
                x=r / r[-1], y=changes[k] + offset, mode="lines",
                name=f"t = {times[k]:.8g} s", line=dict(width=1.5), opacity=0.65,
                hovertemplate="r/R = %{x:.6g}<br>Concentration = %{y:.12g}<extra>%{fullData.name}</extra>"
            ))
profile.add_trace(go.Scatter(
    x=r / r[-1], y=changes[index] + offset, mode="lines",
    name=f"Selected: {times[index]:.8g} s", line=dict(width=3),
    hovertemplate="r/R = %{x:.6g}<br>Concentration = %{y:.12g}<extra>%{fullData.name}</extra>"
))
profile.update_layout(title="Radial concentration profile", height=470,
                      xaxis_title="Fraction of radius, r/R", yaxis_title=axis_title,
                      legend=dict(orientation="h", y=-0.2), margin=dict(l=30, r=25, t=55, b=95))
profile.update_xaxes(range=[0.8 if view["zoom"] else 0, 1])
profile.update_yaxes(tickformat=".10g", exponentformat="none", showexponent="none")
st.plotly_chart(profile, width="stretch", config={"displaylogo": False})

history = go.Figure()
for label, values, dash in [("Surface", changes[:, -1], "solid"),
                             ("Centre", changes[:, 0], "dot"),
                             ("Volume average", averages, "dash")]:
    history.add_trace(go.Scatter(x=times, y=values + offset, mode="lines", name=label,
                                 line=dict(width=2, dash=dash),
                                 hovertemplate="t = %{x:.8g} s<br>Concentration = %{y:.12g}<extra>%{fullData.name}</extra>"))
history.add_vline(x=float(times[index]), line_width=1, line_dash="dot")
history.update_layout(title="Concentration over time", height=340, xaxis_title="Time [s]",
                      yaxis_title=axis_title, legend=dict(orientation="h", y=-0.25),
                      margin=dict(l=30, r=25, t=55, b=80))
history.update_yaxes(tickformat=".10g", exponentformat="none", showexponent="none")
history.update_xaxes(tickformat=".8g")
st.plotly_chart(history, width="stretch", config={"displaylogo": False})

export = io.StringIO()
np.savetxt(export, np.column_stack([r, p["initial"] + changes[index], changes[index]]),
           delimiter=",", header="radius_m,concentration_mol_m3,change_mol_m3", comments="")
st.download_button("Download displayed profile (CSV)", export.getvalue(),
                   file_name=f"concentration_t_{times[index]:.8g}.csv", mime="text/csv")
