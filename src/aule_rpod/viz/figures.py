"""Analysis figures (Part I sensitivity/geometry, attitude characterisation)."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap

from . import style as S

SEQ = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def feasibility_map(ltans, u0s, nwin, path: Path, nominal=(10.5, 0.0), title=None, best=None, D=540.0):
    S.apply()
    fig, ax = plt.subplots(figsize=(11, 4.8))
    from matplotlib.colors import LinearSegmentedColormap
    vmax = max(int(nwin.max()), 1)
    bounds = np.arange(-0.5, vmax + 1.5, 1.0)
    nb = len(bounds) - 1
    if nb <= len(SEQ):
        cmap = ListedColormap(SEQ[:nb])
    else:
        cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ, N=nb)
    norm = BoundaryNorm(bounds, cmap.N)
    dl = (ltans[1] - ltans[0]) / 2
    du = (u0s[1] - u0s[0]) / 2
    im = ax.imshow(nwin, origin="lower", aspect="auto", cmap=cmap, norm=norm,
                   extent=[ltans[0] - dl, ltans[-1] + dl, u0s[0] - du, u0s[-1] + du], interpolation="nearest")
    if best is not None:
        cs = ax.contour(ltans, u0s, best / 60.0, levels=[7.0, 8.0], colors=[S.MUTED, S.INK2], linewidths=0.7)
        ax.clabel(cs, fmt=lambda v: f"best span {v:.0f} min", fontsize=7)
    ax.plot(*nominal, marker="o", ms=9, mfc="none", mec=S.CRITICAL, mew=2)
    ax.annotate("specified orbit\n(LTAN 10:30, u0 = 0 deg)\n0 windows", nominal, xytext=(12, 10), textcoords="offset points",
                fontsize=8, color=S.CRITICAL)
    step = 1 if vmax <= 10 else int(np.ceil(vmax / 10))
    cb = fig.colorbar(im, ax=ax, ticks=range(0, vmax + 1, step), pad=0.01)
    cb.set_label(f"strict docking windows (>= {D / 60:.0f} min) in 10 days")
    cb.outline.set_visible(False)
    ax.set_xlabel("LTAN (local mean solar time of ascending node) [h]")
    ax.set_ylabel("argument of latitude at orbit epoch [deg]")
    ax.grid(False)
    ax.set_title(title or "Window existence vs. the two free orbit parameters", loc="left")
    fig.savefig(path)
    plt.close(fig)
    return path


def relaxation_curves(df, path: Path):
    S.apply()
    studies = [("sun_angle_lower_bound_deg", "illumination lower bound [deg]", True),
               ("uniform_band_widening_deg", "band widening d: [30-d, 60+d] [deg]", False),
               ("min_elevation_deg", "elevation mask [deg]", True),
               ("required_duration_min", "required approach duration [min]", True)]
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.2), sharey=True)
    if "kind" not in df:
        df = df.assign(kind="grid")
    for ax, (key, lab, invert) in zip(axs, studies):
        d = df[(df.study == key) & (df.kind == "grid")].sort_values("value")
        ax.bar(d.value, d.windows, width=0.6 * np.min(np.diff(d.value)) if len(d) > 1 else 0.5, color=S.BLUE)
        for _, m in df[(df.study == key) & (df.kind == "minimum")].iterrows():
            ax.plot([m.value], [m.windows], marker="D", ms=7, color=S.GOOD, mec=S.SURFACE, mew=1.2, zorder=5)
            ax.annotate(f"min {m.value:.2f}: {int(m.windows)}", (m.value, m.windows), xytext=(-7, 0),
                        textcoords="offset points", ha="right", va="center", fontsize=7.5, color=S.INK)
        if not d.windows.any():
            ax.text(0.5, 0.5, "no window at any value", transform=ax.transAxes, ha="center", va="center",
                    fontsize=8.5, color=S.INK2)
        spec = {"sun_angle_lower_bound_deg": 30, "uniform_band_widening_deg": 0, "min_elevation_deg": 5,
                "required_duration_min": 9}[key]
        ax.axvline(spec, color=S.CRITICAL, lw=1.0)
        ax.text(spec, ax.get_ylim()[1] * 0.92 if ax.get_ylim()[1] > 0 else 1, " spec", color=S.CRITICAL, fontsize=8)
        ax.set_xlabel(lab)
        if invert:
            ax.invert_xaxis()
        for x, y in zip(d.value, d.windows):
            if y:
                ax.text(x, y + 0.1, str(int(y)), ha="center", fontsize=7, color=S.INK2)
    axs[0].set_ylabel("strict windows in 10 days")
    fig.suptitle("How much must one requirement give for a window to appear? (others held at spec)",
                 x=0.06, ha="left", fontsize=10.5)
    fig.savefig(path)
    plt.close(fig)
    return path


def attitude_characterisation(sc, law, path: Path, minutes=60.0):
    """Port-axis motion in LVLH (cone about +V-bar) and the sun-angle consequence."""
    S.apply()
    t = sc.t_start + np.arange(0, minutes * 60, 2.0)
    pL = sc.port_axis_lvlh(t)
    fig = plt.figure(figsize=(14, 4.4))
    ax1 = fig.add_subplot(1, 3, 1)
    # gnomonic-like view looking down +Y (V-bar): angular offsets of port from +Y
    off_x = np.degrees(np.arctan2(pL[:, 0], pL[:, 1]))
    off_z = np.degrees(np.arctan2(pL[:, 2], pL[:, 1]))
    ax1.plot(off_z, off_x, color=S.BLUE, lw=0.8)
    ax1.plot([0], [0], marker="+", color=S.INK, ms=10)
    for rr in (5, 10, 15):
        th = np.linspace(0, 2 * np.pi, 200)
        ax1.plot(rr * np.cos(th), rr * np.sin(th), color=S.GRID, lw=0.8)
        ax1.text(rr * 0.71, rr * 0.71, f"{rr} deg", fontsize=7, color=S.MUTED)
    ax1.set_aspect("equal")
    ax1.set_xlabel("toward +Z (orbit normal) [deg]")
    ax1.set_ylabel("toward +X (zenith) [deg]")
    ax1.set_title(f"Docking-port normal seen along +V-bar\n(first {minutes:.0f} min of mission)", loc="left", fontsize=9)
    ax2 = fig.add_subplot(1, 3, (2, 3))
    tt = sc.t_start + np.arange(0, 3 * 3600, 5.0)
    p2 = sc.port_axis_lvlh(tt)
    off = np.degrees(np.arccos(np.clip(p2[:, 1], -1, 1)))
    ax2.plot((tt - sc.t_start) / 60, off, color=S.BLUE, lw=1.0)
    ax2.set_xlabel("minutes from mission start")
    ax2.set_ylabel("port normal offset from +V-bar [deg]")
    ax2.set_ylim(0, 20)
    s = law.summary()
    ax2.set_title("Closed-form law: precession about +Y_LVLH at %.3f deg/s (period %.1f s), nutation %.3f deg,\n"
                  "spin about +X_body at %.6f deg/s (period %.1f s); fit residual %.1e rad over %d samples"
                  % (s["precession_rate_deg_s"], s["precession_period_s"], s["nutation_deg"], s["spin_rate_deg_s"],
                     s["spin_period_s"], s["max_residual_rad"], len(sc.profile.times_s)), loc="left", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def geometry_explainer(sc, cs, path: Path):
    """Why only Wellington can host a window: sun-angle band vs latitude on the day arc,
    with each station's coverage on a separate panel (no mixed units)."""
    S.apply()
    T = sc.orbit.nodal_period_s if hasattr(sc.orbit, "nodal_period_s") else 5746.0
    t = sc.t_start + np.arange(0, 3 * T, 5.0)
    r, v = sc.state_gcrs(t)
    from ..frames import itrs_to_geodetic
    lat, lon, _ = itrs_to_geodetic(np.einsum("nij,nj->ni", sc.eo.gcrs_to_itrs(t), r), sc.re, sc.f)
    sL = sc.sun_dir_lvlh(t)
    th_y = np.degrees(np.arccos(np.clip(sL[:, 1], -1, 1)))
    ang = sc.sun_angle_deg(t)
    lit = sc.g_full_sun(t) >= 0
    asc = np.gradient(lat) > 0
    sel = asc & lit
    fig, (ax, axb) = plt.subplots(2, 1, figsize=(11, 6.2), sharex=True, gridspec_kw={"height_ratios": [3.0, 1.15]})
    ax.axhspan(cs.sun_min_deg, cs.sun_max_deg, color=S.BAND_COLOR, alpha=0.15, lw=0)
    ax.scatter(lat[sel], ang[sel], s=1.5, color=S.MUTED, label="Sun angle to port normal (3 orbits)")
    ax.scatter(lat[sel], th_y[sel], s=1.5, color=S.INK, label="Sun angle to +V-bar = cone axis of the port")
    ax.set_ylim(0, 180)
    ax.set_ylabel("angle [deg]")
    ax.legend(loc="upper left", markerscale=6)
    ax.set_title("Ascending (daylight) arc: the 30-60 deg band holds only south of ~12 S, where only Wellington looks",
                 loc="left", fontsize=10)
    lam = np.degrees(np.arccos(sc.re * np.cos(np.radians(5)) / np.linalg.norm(r[0])) - np.radians(5))
    core = 7.4    # central-angle radius that still gives a >= 9-min pass (550 km, 5 deg mask)
    names = [s_.name for s_ in sc.stations]
    for k, st in enumerate(sc.stations):
        y = len(names) - 1 - k
        axb.plot([st.lat_deg - lam, min(st.lat_deg + lam, 90)], [y, y], color=S.STATION_COLORS[st.name], lw=7,
                 solid_capstyle="butt")
        axb.plot([st.lat_deg - core, min(st.lat_deg + core, 90)], [y, y], color=S.INK, lw=1.6)
    axb.set_yticks(range(len(names)))
    axb.set_yticklabels(names[::-1], fontsize=8)
    axb.set_ylim(-0.7, len(names) - 0.3)
    axb.grid(axis="y", visible=False)
    axb.set_xlim(-90, 90)
    axb.set_xlabel("target latitude [deg]   (coloured: visible at >= 5 deg elevation; black: overhead enough for a >= 9-min pass)")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path
