"""Gantt charts of constraint satisfaction and docking windows (static + interactive)."""
from __future__ import annotations

from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..events import intervals as iv
from . import style as S


def _to_dt(sc, t):
    return pd.to_datetime(sc.ts.datetime(t))


def _rows(sc, cs, duration_s):
    rows = []
    for st in sc.stations:
        P = cs.passes[st.name]
        long_ = iv.filter_min_duration(P, duration_s)
        short = np.array([p for p in P if p[1] - p[0] < duration_s - 1e-9]).reshape(-1, 2)
        rows.append((f"{st.name}\nGS pass (el >= {cs.min_el_deg:g} deg)", [(long_, S.STATION_COLORS.get(st.name, S.BLUE), 1.0),
                                                                     (short, S.STATION_COLORS.get(st.name, S.BLUE), 0.30)]))
    rows.append(("Full sunlight\n(no umbra/penumbra)", [(cs.lit, S.SUN_COLOR, 0.85)]))
    rows.append((f"Sun angle in\n[{cs.sun_min_deg:g}, {cs.sun_max_deg:g}] deg", [(cs.band, S.BAND_COLOR, 0.9)]))
    lb = iv.intersect(cs.lit, cs.band)
    rows.append(("Sunlight AND band\n(>= 9 min shown dark)", [(iv.filter_min_duration(lb, duration_s), S.INK2, 1.0),
                                                             (lb, S.MUTED, 0.35)]))
    return rows


def gantt_overview(sc, cs, windows, duration_s, path: Path, relaxed=None, candidates=None, title=None):
    S.apply()
    rows = _rows(sc, cs, duration_s)
    n = len(rows) + 1
    fig, ax = plt.subplots(figsize=(16, 0.62 * n + 1.6))
    day = 86400.0
    t0 = sc.t_start
    for k, (label, layers) in enumerate(rows):
        y = n - 1 - k
        for ivs, color, alpha in layers:
            if len(ivs):
                ax.broken_barh([((a - t0) / day, max((b - a) / day, 1e-4)) for a, b in ivs], (y - 0.32, 0.64),
                               facecolors=color, edgecolors=color, linewidth=0.35, alpha=alpha)
    # docking windows row
    y = 0
    if len(windows):
        ax.broken_barh([((w.feasible_start_s - t0) / day, (w.feasible_end_s - w.feasible_start_s) / day)
                        for w in windows], (y - 0.32, 0.64), facecolors=S.GOOD, edgecolors=S.GOOD, linewidth=0.6)
    else:
        ax.text(0.05, y + 0.52, "X  NO STRICT WINDOW: 0 of %d passes >= 9 min satisfy all constraints for %.0f min "
                "(orange: best feasible span of the closest passes, ranked)"
                % (sum(len(iv.filter_min_duration(P, duration_s)) for P in cs.passes.values()), duration_s / 60),
                ha="left", va="center", fontsize=8, color=S.CRITICAL, fontweight="bold")
    if candidates is not None and len(candidates):
        for _, c in candidates.head(6).iterrows():
            ts_ = sc.ts.t_from_iso(c["longest_feasible_start_utc"].rstrip("Z"))
            te_ = sc.ts.t_from_iso(c["longest_feasible_end_utc"].rstrip("Z"))
            xm = (ts_ - t0) / day
            ax.plot([xm, (te_ - t0) / day], [y, y], color=S.SERIOUS, lw=9, solid_capstyle="butt")
            ax.annotate(f"#{int(c['rank'])} {c['longest_feasible_s'] / 60:.1f} min", (xm, y), xytext=(0, 9),
                        textcoords="offset points", ha="center", fontsize=7, color=S.INK2)
    if relaxed:
        for lab, ws in relaxed.items():
            for w in ws:
                xm = (w.feasible_start_s - t0) / day
                ax.plot([xm], [y - 0.05], marker="D", ms=7, color=S.GOOD, mec=S.SURFACE, mew=1.2, zorder=5)
                ax.annotate(lab, (xm, y), xytext=(0, -17), textcoords="offset points", ha="center",
                            fontsize=7, color=S.INK2)
    labels = [r[0] for r in rows] + ["Docking window\n(all constraints >= 9 min)"]
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels[::-1], fontsize=8, color=S.INK)
    ax.set_ylim(-0.9, n - 0.4)
    ax.set_xlim(0, (sc.t_end - t0) / day)
    ticks = np.arange(0, (sc.t_end - t0) / day + 1e-9, 1.0)
    ax.set_xticks(ticks)
    ax.set_xticklabels([sc.ts.iso(t0 + d * day, 0)[0][5:16].replace("T", "\n") for d in ticks], fontsize=7.5)
    ax.set_xticks(np.arange(0, ticks[-1], 0.25), minor=True)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", which="minor", color=S.GRID, linewidth=0.3)
    ax.set_xlabel("UTC (mission window %s to %s)" % (sc.ts.iso(t0, 0)[0], sc.ts.iso(sc.t_end, 0)[0]))
    ax.set_title(title or "Docking-window constraint Gantt - 10-day mission window", loc="left", fontsize=11, color=S.INK)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Patch(color=S.STATION_COLORS[s.name], label=f"{s.name} pass >= 9 min") for s in sc.stations if s.name in S.STATION_COLORS]
    handles += [Patch(color=S.MUTED, alpha=0.35, label="pass < 9 min (faint)"),
                Patch(color=S.SUN_COLOR, label="full sunlight"), Patch(color=S.BAND_COLOR, label="sun angle in band"),
                Line2D([], [], color=S.SERIOUS, lw=4, label="best near-miss (longest feasible < 9 min)"),
                Line2D([], [], marker="D", color=S.GOOD, lw=0, label="window under minimum relaxation")]
    ax.legend(handles=handles, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=7.5)
    fig.savefig(path)
    fig.savefig(Path(path).with_suffix(".pdf"))
    plt.close(fig)
    return path


def candidate_panels(sc, cs, candidates, duration_s, path: Path, n=6, pad_s=90.0):
    """Anatomy of the closest passes: sun-angle trace + mini-Gantt per pass."""
    S.apply()
    cands = candidates.head(n)
    ncol = 2
    nrow = int(np.ceil(len(cands) / ncol))
    fig = plt.figure(figsize=(16, 4.6 * nrow))
    outer = fig.add_gridspec(nrow, ncol, hspace=0.42, wspace=0.22, top=0.95, bottom=0.05, left=0.07, right=0.98)
    for idx, (_, c) in enumerate(cands.iterrows()):
        r, col = divmod(idx, ncol)
        inner = outer[r, col].subgridspec(2, 1, height_ratios=[1.6, 1.15], hspace=0.08)
        aos = sc.ts.t_from_iso(c["pass_aos_utc"].rstrip("Z"))
        los = sc.ts.t_from_iso(c["pass_los_utc"].rstrip("Z"))
        ta, tb = aos - pad_s, aos + 13.0 * 60.0
        t = np.arange(ta, tb, 1.0)
        ang = sc.sun_angle_deg(t)
        sun_l = sc.sun_dir_lvlh(t)
        cone_axis = np.degrees(np.arccos(np.clip(sun_l[:, 1], -1, 1)))
        x = (t - aos) / 60.0
        ax1 = fig.add_subplot(inner[0])
        ax2 = fig.add_subplot(inner[1], sharex=ax1)
        ax1.axhspan(cs.sun_min_deg, cs.sun_max_deg, color=S.BAND_COLOR, alpha=0.13, lw=0)
        ax1.plot(x, cone_axis, color=S.MUTED, lw=1.0, label="Sun angle to +V-bar (port cone axis)")
        ax1.plot(x, ang, color=S.INK, lw=1.4, label="Sun angle to docking-port normal")
        for yb in (cs.sun_min_deg, cs.sun_max_deg):
            ax1.axhline(yb, color=S.BAND_COLOR, lw=0.8)
        ax1.set_ylabel("sun angle [deg]")
        ax1.set_ylim(0, 90)
        if idx == 0:
            ax1.legend(loc="upper left", fontsize=7)
        ax1.set_title(f"#{int(c['rank'])}  {c['station']}  AOS {c['pass_aos_utc'][:19].replace('T', ' ')}Z   "
                      f"max el {c['max_elevation_deg']:.1f} deg   longest feasible {c['longest_feasible_s'] / 60:.2f} min "
                      f"(need {duration_s / 60:.0f})", loc="left", fontsize=8.5)
        plt.setp(ax1.get_xticklabels(), visible=False)
        rows = [("GS pass", [[aos, los]], S.STATION_COLORS.get(c["station"], S.BLUE)),
                ("Full sunlight", iv.intersect(cs.lit, [[ta, tb]]), S.SUN_COLOR),
                (f"sun >= {cs.sun_min_deg:g} deg", iv.intersect(cs.sun_ge_min, [[ta, tb]]), S.BAND_COLOR),
                (f"sun <= {cs.sun_max_deg:g} deg", iv.intersect(cs.sun_le_max, [[ta, tb]]), S.BAND_COLOR),
                ("ALL constraints", iv.intersect_all([[aos, los]], cs.lit, cs.band), S.SERIOUS)]
        for k, (lab, ivs, colr) in enumerate(rows):
            yy = len(rows) - 1 - k
            if len(ivs):
                ax2.broken_barh([((a - aos) / 60, (b - a) / 60) for a, b in ivs], (yy - 0.3, 0.6), facecolors=colr, alpha=0.9)
        ls = sc.ts.t_from_iso(c["longest_feasible_start_utc"].rstrip("Z")) if c["longest_feasible_start_utc"] else aos
        req = ((ls - aos) / 60, (ls - aos + duration_s) / 60)
        ax2.annotate("", xy=(req[1], -0.8), xytext=(req[0], -0.8), annotation_clip=False,
                     arrowprops=dict(arrowstyle="<->", color=S.CRITICAL, lw=1.0))
        ax2.text(np.mean(req), -1.55, "required %.0f-min approach overruns the feasible span by %.0f s"
                 % (duration_s / 60, c["deficit_vs_required_s"]), ha="center", fontsize=7, color=S.CRITICAL)
        ax2.axvline((los - aos) / 60, color=S.AXIS, lw=0.8)
        ax2.set_yticks(range(len(rows)))
        ax2.set_yticklabels([rw[0] for rw in rows][::-1], fontsize=7)
        ax2.set_ylim(-2.0, len(rows) - 0.4)
        ax2.grid(axis="y", visible=False)
        ax2.set_xlabel("minutes from AOS", labelpad=1)
        ax2.set_xlim(x[0], x[-1])
    fig.suptitle("Closest docking opportunities: the port's 5-min nutation rides on the orbital ramp, "
                 "so no 9-min span stays inside the 30-60 deg band", x=0.07, ha="left", fontsize=11, y=0.995)
    fig.savefig(path)
    plt.close(fig)
    return path


def gantt_interactive(sc, cs, windows, duration_s, path: Path, relaxed=None, candidates=None):
    import plotly.express as px
    import plotly.graph_objects as go
    recs = []

    def add(label, ivs, cat, detail=""):
        for a, b in ivs:
            recs.append({"Constraint": label, "Start": _to_dt(sc, a)[0], "Finish": _to_dt(sc, b)[0],
                         "Category": cat, "Duration (min)": round((b - a) / 60, 2), "Detail": detail})
    for st in sc.stations:
        P = cs.passes[st.name]
        add(f"{st.name} pass", iv.filter_min_duration(P, duration_s), f"{st.name} pass >= 9 min")
        add(f"{st.name} pass", np.array([p for p in P if p[1] - p[0] < duration_s]).reshape(-1, 2), "pass < 9 min")
    add("Full sunlight", cs.lit, "full sunlight")
    add("Sun angle 30-60 deg", cs.band, "sun angle in band")
    lb = iv.intersect(cs.lit, cs.band)
    add("Sunlight AND band", lb, "sunlight AND band")
    add("Docking window (strict)", [(w.feasible_start_s, w.feasible_end_s) for w in windows], "DOCKING WINDOW")
    if candidates is not None:
        for _, c in candidates.iterrows():
            if not c["longest_feasible_start_utc"]:
                continue
            a = sc.ts.t_from_iso(c["longest_feasible_start_utc"].rstrip("Z"))
            b = sc.ts.t_from_iso(c["longest_feasible_end_utc"].rstrip("Z"))
            add("Best feasible span per pass", [(a, b)], "near-miss (< 9 min)",
                f"rank {int(c['rank'])}; deficit {c['deficit_vs_required_s']:.0f} s; limited by {c['limited_at_start_by']} / {c['limited_at_end_by']}")
    if relaxed:
        for lab, ws in relaxed.items():
            add(f"Relaxed: {lab}", [(w.feasible_start_s, w.feasible_end_s) for w in ws], "relaxed window")
    df = pd.DataFrame(recs)
    cmap = {f"{s.name} pass >= 9 min": S.STATION_COLORS.get(s.name, S.BLUE) for s in sc.stations}
    cmap.update({"pass < 9 min": "#c3c2b7", "full sunlight": S.SUN_COLOR, "sun angle in band": S.BAND_COLOR,
                 "sunlight AND band": S.INK2, "DOCKING WINDOW": S.GOOD, "near-miss (< 9 min)": S.SERIOUS,
                 "relaxed window": S.GREEN})
    order = [f"{s.name} pass" for s in sc.stations] + ["Full sunlight", "Sun angle 30-60 deg", "Sunlight AND band",
                                                         "Best feasible span per pass"] + \
            [f"Relaxed: {k}" for k in (relaxed or {})] + ["Docking window (strict)"]
    fig = px.timeline(df, x_start="Start", x_end="Finish", y="Constraint", color="Category",
                      color_discrete_map=cmap, hover_data=["Duration (min)", "Detail"],
                      category_orders={"Constraint": order})
    fig.update_traces(marker_line_width=0.6, marker_line_color="rgba(0,0,0,0.25)")
    if not windows:
        fig.add_annotation(x=_to_dt(sc, sc.t_start + 5 * 86400)[0], y="Docking window (strict)",
                           text="NO STRICT DOCKING WINDOW (0 passes satisfy all constraints for 9 min)",
                           showarrow=False, font=dict(color=S.CRITICAL, size=12))
    fig.update_layout(title="Docking-window constraint Gantt (zoom: drag; reset: double-click)",
                      plot_bgcolor=S.SURFACE, paper_bgcolor=S.SURFACE, height=640, legend_title_text="",
                      xaxis=dict(rangeslider=dict(visible=True), type="date"),
                      font=dict(family="system-ui, -apple-system, Segoe UI, sans-serif", size=12))
    fig.update_yaxes(autorange="reversed")
    fig.write_html(str(path), include_plotlyjs=True, full_html=True)
    return path
