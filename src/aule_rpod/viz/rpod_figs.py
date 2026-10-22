"""Part II figures: approach trajectory, timeline, strategy trade, safety, dispersions."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from . import style as S

PHASE_COLORS = {"leg 1 transfer": S.BLUE, "leg 2 corridor": S.VIOLET, "hold (go/no-go)": S.ORANGE,
                "final translation": S.GREEN, "final braking/contact": S.GREEN, "contact": S.GREEN}


def _phase_lines(ax, plan, xi, yi, lw=2.0, label=True):
    for ph, col in PHASE_COLORS.items():
        k = plan.phase_of_t == ph
        if k.any() and ph != "contact":
            ax.plot(plan.r[k, xi], plan.r[k, yi], color=col, lw=lw, label=ph if label else None)


def trajectory(plan, path: Path, kos_m, corridor_deg, sun_lvlh_at_hp=None, drift=None, title=None, r_col=5.0):
    S.apply()
    fig = plt.figure(figsize=(16, 9.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1.25, 1.0], hspace=0.30, wspace=0.16)
    th = np.linspace(0, 2 * np.pi, 300)
    # ---------------- (a) full in-plane view
    ax = fig.add_subplot(gs[0, 0])
    ax.fill(kos_m * np.cos(th), kos_m * np.sin(th), color=S.GRID, alpha=0.6, lw=0)
    ax.plot(kos_m * np.cos(th), kos_m * np.sin(th), color=S.AXIS, lw=0.8)
    ax.text(0, -kos_m - 5, f"keep-out sphere R = {kos_m:g} m\n(corridor-only inside)", fontsize=7.5, color=S.INK2,
            ha="center", va="top")
    if drift:
        for d in drift:
            X = d["trajectory"]
            k = (np.abs(X[:, 1]) < 140) & (X[:, 0] > -95) & (X[:, 0] < 40)
            ax.plot(X[k, 1], X[k, 0], color=S.MUTED, lw=0.9)
    _phase_lines(ax, plan, 1, 0)
    b1 = plan.burns[0]
    ax.plot(b1.r[1], b1.r[0], "o", ms=7, color=S.INK, mec=S.SURFACE, mew=1.2, zorder=5)
    ax.annotate(f"P0 / B1  {b1.mag * 100:.1f} cm/s\n(t = 0, window start)", (b1.r[1], b1.r[0]), xytext=(-150, -6),
                textcoords="offset points", fontsize=8)
    ax.plot(0, 0, "s", ms=8, color=S.INK)
    if sun_lvlh_at_hp is not None:
        s_ = sun_lvlh_at_hp / np.linalg.norm(sun_lvlh_at_hp[:2])
        ax.annotate("", xy=(30 + 35 * s_[1], 10 + 35 * s_[0]), xytext=(30, 10),
                    arrowprops=dict(arrowstyle="->", color=S.YELLOW, lw=2.2))
        ax.text(30 + 36 * s_[1], 10 + 36 * s_[0], "Sun (in-plane projection,\nduring hold)", color=S.INK2, fontsize=7.5)
    ax.annotate("", xy=(136, 16), xytext=(136, 2), arrowprops=dict(arrowstyle="->", color=S.MUTED, lw=1.0))
    ax.text(134, 9, "zenith (+R)", fontsize=7, color=S.MUTED, ha="right", va="center")
    ax.annotate("", xy=(138, -4), xytext=(116, -4), arrowprops=dict(arrowstyle="->", color=S.MUTED, lw=1.0))
    ax.text(114, -4, "orbital velocity (+V) ", fontsize=7, color=S.MUTED, ha="right", va="center")
    if drift:
        ax.text(-22, 30, "grey: free drift if a burn is missed\n(B2 missed: min range %.1f m)"
                % [d for d in drift if d["failed_burn"] == "B2"][0]["min_range_m"], fontsize=7.5, color=S.INK2, va="top")
    ax.set_xlabel("y  along-track / V-bar [m]   (+ = ahead of target)")
    ax.set_ylabel("x  radial / R-bar [m]   (+ = up)")
    ax.set_aspect("equal")
    ax.set_xlim(-25, 140)
    ax.set_ylim(-95, 38)
    ax.legend(loc="lower left", fontsize=7.5)
    ax.set_title(title or "Approach in the orbital plane (target LVLH)", loc="left")
    # ---------------- (b) near-field zoom
    az = fig.add_subplot(gs[0, 1])
    L = kos_m * 1.08
    ca, sa = np.cos(np.radians(corridor_deg)), np.sin(np.radians(corridor_deg))
    az.fill([0, L * ca, L * ca], [0, L * sa, -L * sa], color=S.GREEN, alpha=0.10, lw=0)
    az.plot([0, L * ca], [0, L * sa], color=S.GREEN, lw=0.7)
    az.plot([0, L * ca], [0, -L * sa], color=S.GREEN, lw=0.7)
    az.plot(kos_m * np.cos(th), kos_m * np.sin(th), color=S.AXIS, lw=0.8)
    az.plot(r_col * np.cos(th), r_col * np.sin(th), color=S.CRITICAL, lw=0.8)
    az.text(-r_col * 0.2, -r_col - 1.2, f"collision envelope {r_col:g} m", fontsize=7, color=S.CRITICAL, ha="center", va="top")
    if drift:
        d3 = [d for d in drift if d["failed_burn"] == "B3"]
        for d in d3:
            X = d["trajectory"]
            k = (np.abs(X[:, 1]) < 30) & (np.abs(X[:, 0]) < 16)
            az.plot(X[k, 1], X[k, 0], color=S.MUTED, lw=0.9)
            az.text(-14, -12.5, f"grey: B3 missed -> min range {d['min_range_m']:.1f} m\n(covered by the posigrade\nbreakout burn)",
                    fontsize=7, color=S.INK2)
    _phase_lines(az, plan, 1, 0, label=False)
    lab = {"B2": (6, 10), "B3": (-6, 14), "B4": (-70, 16)}
    for b in plan.burns[1:]:
        az.plot(b.r[1], b.r[0], "o", ms=6, color=S.INK, mec=S.SURFACE, mew=1.2, zorder=5)
        az.annotate(f"{b.id}  {b.mag * 100:.1f} cm/s", (b.r[1], b.r[0]), xytext=lab.get(b.id, (5, 5)),
                    textcoords="offset points", fontsize=8,
                    arrowprops=dict(arrowstyle="-", color=S.MUTED, lw=0.6))
    az.text(plan.params["HP"][1] + 0.5, -3.6, "10 m hold\n(go/no-go)", fontsize=7.5, ha="center", color=S.ORANGE)
    az.text(kos_m - 1, 5.2, "corridor entry\n(WP1, KOS boundary)", fontsize=7.5, ha="center", color=S.VIOLET)
    az.text(L * 0.55, -L * sa * 0.55 - 2.6, f"approach corridor +/-{corridor_deg:g} deg", fontsize=7.5, color=S.GREEN)
    az.plot(0, 0, "s", ms=9, color=S.INK)
    az.text(-0.8, 1.2, "port", fontsize=8, ha="right")
    az.set_aspect("equal")
    az.set_xlim(-16, 31)
    az.set_ylim(-16, 16)
    az.set_xlabel("y along-track [m]")
    az.set_ylabel("x radial [m]")
    az.set_title("Near field: corridor, hold point and burns", loc="left", fontsize=9.5)
    # ---------------- (c) cross-track
    ax2 = fig.add_subplot(gs[1, 0])
    _phase_lines(ax2, plan, 1, 2, label=False)
    ax2.plot(0, 0, "s", ms=7, color=S.INK)
    ax2.set_xlabel("y along-track [m]")
    ax2.set_ylabel("z cross-track [m]")
    ax2.set_title("Cross-track view: the 20 m out-of-plane offset is nulled on leg 1", loc="left", fontsize=9.5)
    # ---------------- (d) range & closing speed
    ax3 = fig.add_subplot(gs[1, 1])
    rng = np.linalg.norm(plan.r, axis=1)
    for ph, col in PHASE_COLORS.items():
        k = plan.phase_of_t == ph
        if k.any() and ph != "contact":
            ax3.plot(plan.t[k], rng[k], color=col, lw=1.8)
    ax3.set_yscale("log")
    ax3.set_ylim(0.5, 200)
    ax3.set_xlabel("time from window start [s]")
    ax3.set_ylabel("range to port [m] (log)")
    cs_ = -np.sum(plan.r * plan.v, 1) / np.maximum(rng, 1e-9)
    leg1 = plan.phase_of_t == "leg 1 transfer"
    ax3.set_title("Range history: closing %.2f m/s (leg 1) -> %.2f m/s (contact)"
                  % (np.nanmax(cs_[leg1]), plan.params["v_contact"]), loc="left", fontsize=9.5)
    fig.savefig(path)
    plt.close(fig)
    return path


def timeline(plan, path: Path, t0_abs, window_s, sun_on_port_fn, gs_name, band=(30.0, 60.0), utc_fn=None):
    S.apply()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 5.6), sharex=True, gridspec_kw={"height_ratios": [1.2, 1.0]})
    rows = [(f"{gs_name} pass", [(0, window_s)], S.AQUA), ("Full sunlight", [(0, window_s)], S.YELLOW)]
    for ph, t_a, t_b in plan.phases:
        rows.append((ph, [(t_a, t_b)], PHASE_COLORS.get(ph, S.MUTED)))
    for k, (lab, ivs, col) in enumerate(rows):
        y = len(rows) - 1 - k
        ax1.broken_barh([(a, b - a) for a, b in ivs], (y - 0.3, 0.6), facecolors=col)
        ax1.text(-4, y, lab, ha="right", va="center", fontsize=8)
    for b in plan.burns:
        ax1.axvline(b.t_s, color=S.INK, lw=0.6)
        ax1.text(b.t_s, len(rows) - 0.35, b.id, fontsize=7.5, ha="center", va="bottom")
    ax1.axvline(plan.t_contact, color=S.GOOD, lw=1.2)
    ax1.text(plan.t_contact - 4, -0.42, "contact %.0f s\n(margin %.0f s to LOS)" % (plan.t_contact, window_s - plan.t_contact),
             fontsize=7.5, color=S.GOOD, ha="right", va="top")
    ax1.set_yticks([])
    ax1.set_ylim(-1.3, len(rows) + 0.2)
    ax1.grid(axis="y", visible=False)
    ttl = "Part II timeline inside the docking window"
    if utc_fn:
        ttl += f"  (t = 0 at {utc_fn(t0_abs)})"
    ax1.set_title(ttl, loc="left")
    tt = np.linspace(0, window_s, 400)
    ang = sun_on_port_fn(t0_abs + tt)
    ax2.axhspan(*band, color=S.MAGENTA, alpha=0.13, lw=0)
    ax2.plot(tt, ang, color=S.INK, lw=1.6, label="Sun angle on docking port (port held on +V-bar)")
    t_hp = plan.params["T1"] + plan.params["T2"]
    ax2.axvspan(t_hp, plan.t_contact, color=S.GREEN, alpha=0.08, lw=0)
    ax2.text((t_hp + plan.t_contact) / 2, band[1] + 3, "critical phase (hold + final translation)", ha="center",
             fontsize=8, color=S.GREEN)
    ax2.set_ylabel("Sun angle [deg]")
    ax2.set_xlabel("time from window start [s]")
    ax2.set_ylim(0, 90)
    ax2.legend(loc="lower right")
    ax2.set_xlim(-5, window_s + 5)
    fig.subplots_adjust(left=0.13)
    fig.savefig(path)
    plt.close(fig)
    return path


def trade_small_multiples(df, plans, sun_fn, t0_abs, path: Path, band=(30.0, 60.0)):
    """Sun angle on the re-oriented port for every candidate axis: selected axis first,
    then the other eligible ones by rank, then the rejected ones by delta-v."""
    S.apply()
    d = df.copy()
    d["_el"] = d.get("eligible", False).astype(bool)
    d["_rank"] = d.get("rank").fillna(99) if "rank" in d else 99
    d = d.sort_values(["_el", "_rank", "dv_total_mps"], ascending=[False, True, True])
    keys = [k for k in d.strategy if k in plans]
    ncol = 4
    nrow = int(np.ceil((len(keys) + 1) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(12.5, 3.1 * nrow), sharey=True, sharex=True)
    axs = np.atleast_1d(axs).ravel()
    for ax, k in zip(axs, keys):
        p = plans[k]
        row = d[d.strategy == k].iloc[0]
        tt = np.linspace(0, p.t_contact, 200)
        ang = np.degrees(np.arccos(np.clip(sun_fn(t0_abs + tt) @ p.strategy.axis, -1, 1)))
        ax.axhspan(*band, color=S.MAGENTA, alpha=0.13, lw=0)
        t_hp = p.params["T1"] + p.params["T2"]
        ax.axvspan(t_hp, p.t_contact, color=S.GREEN, alpha=0.08, lw=0)
        ax.plot(tt, ang, color=S.INK, lw=1.5)
        el = bool(row["_el"])
        verdict = ("selected" if row["_rank"] == 1 else "eligible") if el else "rejected"
        ax.set_title(f"{row['name']}\n\u0394v {row['dv_total_mps']:.2f} m/s, {verdict}",
                     fontsize=8.5, loc="left", color=S.INK if el else S.CRITICAL)
        ax.set_ylim(0, 180)
        ax.set_yticks([0, 30, 60, 90, 120, 150, 180])
    for ax in axs[len(keys):]:
        ax.axis("off")
    if len(axs) > len(keys):
        axs[len(keys)].text(0.02, 0.62, "shaded band: 30-60 deg requirement\ngreen: hold + final translation\n"
                            "red title: rejected (dark port, out of band\nor Earth in the camera background)",
                            transform=axs[len(keys)].transAxes, fontsize=8.5, color=S.INK2, va="top")
    for ax in axs[::ncol]:
        ax.set_ylabel("Sun angle on port [deg]")
    for ax in axs[(nrow - 1) * ncol:]:
        ax.set_xlabel("time from window start [s]")
    for ax in axs[len(keys) - ncol:len(keys)]:
        ax.xaxis.set_tick_params(labelbottom=True)
        ax.set_xlabel("time from window start [s]")
    fig.suptitle("Candidate approach axes at the selected window: Sun angle on the re-oriented docking port",
                 x=0.01, ha="left", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def dispersions(mc_raw, plan, path: Path, corridor_deg):
    """Closed-loop Monte Carlo: hold-point dispersion, the three translational IDSS
    initial-contact quantities and the delta-v distribution."""
    from ..rpod.montecarlo import IDSS_CLOSING_RATE_MPS, IDSS_LATERAL_MISALIGNMENT_M, IDSS_LATERAL_RATE_MPS
    S.apply()
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 6.6))

    def hist(ax, x, color, xlabel, p99_fmt=None, limit=None, band=None, nominal=None):
        x = np.asarray(x, float)
        x = x[np.isfinite(x)]
        if band is not None:
            ax.axvspan(*band, color=S.BAND_COLOR, alpha=0.18, lw=0)
        ax.hist(x, bins=60, color=color, edgecolor=S.SURFACE, linewidth=0.3)
        top = ax.get_ylim()[1]
        if p99_fmt:
            q = np.percentile(x, 99)
            ax.axvline(q, color=S.INK, lw=1.0)
            ax.text(q, top * 0.92, " p99 " + p99_fmt % q, color=S.INK, fontsize=8, va="top")
        if limit is not None:
            ax.axvline(limit[0], color=S.CRITICAL, lw=1.2, ls="--")
            ax.text(limit[0], top * 0.72, limit[1], color=S.CRITICAL, fontsize=8, ha="right", va="top")
        if nominal is not None:
            ax.axvline(nominal, color=S.MUTED, lw=1.0, ls=":")
            ax.text(nominal, top * 0.92, " nominal %.3f m/s" % nominal, color=S.MUTED, fontsize=8, ha="left", va="top")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("runs")

    hist(axs[0, 0], mc_raw["hp_pos_err"], S.BLUE, "arrival error at the 10-m hold point [m]", "%.2f m")
    hist(axs[0, 1], np.asarray(mc_raw["lateral_at_contact"]) * 100, S.AQUA, "lateral misalignment at contact [cm]", "%.1f cm",
         limit=(IDSS_LATERAL_MISALIGNMENT_M * 100, "IDSS limit 10 cm "))
    xc = np.asarray(mc_raw["contact_speed"]) * 100
    hist(axs[1, 0], xc, S.VIOLET, "closing speed at contact [cm/s]",
         limit=(IDSS_CLOSING_RATE_MPS[0] * 100, "IDSS min 5 cm/s "))
    hi = max(7.0, np.nanmax(xc) + 0.3)
    axs[1, 0].set_xlim(4.6, hi)
    top = axs[1, 0].get_ylim()[1]
    axs[1, 0].text(hi, top * 0.92, "IDSS max %.0f cm/s \u2192" % (IDSS_CLOSING_RATE_MPS[1] * 100),
                   color=S.CRITICAL, fontsize=8, ha="right", va="top")
    axs[1, 0].text(hi, top * 0.80, "p1-p99  %.2f-%.2f cm/s" % tuple(np.nanpercentile(xc, [1, 99])),
                   color=S.INK, fontsize=8, ha="right", va="top")
    hist(axs[1, 1], mc_raw["dv_total"], S.ORANGE, "total delta-v [m/s]", "%.3f m/s", nominal=plan.dv_total)
    n = len(np.asarray(mc_raw["dv_total"]))
    lr = 1000 * np.nanpercentile(np.asarray(mc_raw["contact_lateral_rate"], float), 99)
    ok = 100 * np.mean(np.asarray(mc_raw["idss_ok"], bool))
    fig.suptitle(f"Closed-loop Monte Carlo, {n:,} runs. Contact: {ok:.1f} % inside the IDSS translational envelope; "
                 f"lateral rate p99 {lr:.1f} mm/s (limit {IDSS_LATERAL_RATE_MPS * 1000:.0f} mm/s)", x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def timeline_trade(rows, chosen, path: Path, b3_max=0.15, miss_min=10.0):
    """Split of the transfer time between the two legs: total delta-v falls with a longer
    first leg while the stop burn B3 grows and the miss distance if B2 fails shrinks.
    The chosen split is the cheapest one that meets both limits."""
    S.apply()
    T1 = np.array([r["T1"] for r in rows])
    dv = np.array([r["dv_total"] for r in rows])
    b3 = np.array([r["B3"] for r in rows]) * 100
    miss = np.array([r["miss_B2"] for r in rows])
    ok = (b3 <= b3_max * 100) & (miss >= miss_min)
    fig, axs = plt.subplots(1, 3, figsize=(14, 3.7), sharex=True)
    for ax in axs:
        if ok.any():
            ax.axvspan(T1[ok].min(), T1[ok].max(), color=S.GREEN, alpha=0.07, lw=0)
        ax.axvline(chosen["T1"], color=S.CRITICAL, lw=0.8, ls=":")
        ax.set_xlabel("leg-1 duration T1 [s]  (T2 = budget - T1)")
    axs[0].plot(T1, dv, color=S.BLUE, lw=1.7)
    axs[0].plot(chosen["T1"], chosen["dv_total"], "o", color=S.CRITICAL)
    axs[0].annotate(" chosen: T1 = %.0f s, %.3f m/s" % (chosen["T1"], chosen["dv_total"]), (chosen["T1"], chosen["dv_total"]),
                    xytext=(6, 8), textcoords="offset points", fontsize=8, color=S.INK)
    axs[0].set_ylabel("total delta-v [m/s]")
    axs[0].set_title("Total delta-v (minimised)", loc="left", fontsize=9.5)
    axs[1].plot(T1, b3, color=S.ORANGE, lw=1.7)
    axs[1].axhline(b3_max * 100, color=S.CRITICAL, lw=1.0, ls="--")
    axs[1].text(T1.min(), b3_max * 100, " limit %.0f cm/s" % (b3_max * 100), color=S.CRITICAL, fontsize=8, va="bottom")
    axs[1].plot(chosen["T1"], chosen["B3"] * 100, "o", color=S.CRITICAL)
    axs[1].set_ylabel("stop burn at the hold point, B3 [cm/s]")
    axs[1].set_title("Braking impulse at 10 m (plume exposure)", loc="left", fontsize=9.5)
    axs[2].plot(T1, miss, color=S.VIOLET, lw=1.7)
    axs[2].axhline(miss_min, color=S.CRITICAL, lw=1.0, ls="--")
    axs[2].text(T1.max(), miss_min, "limit %.0f m " % miss_min, color=S.CRITICAL, fontsize=8, va="bottom", ha="right")
    axs[2].plot(chosen["T1"], chosen["miss_B2"], "o", color=S.CRITICAL)
    axs[2].set_ylabel("closest approach if B2 is missed [m]")
    axs[2].set_title("Passive safety of the corridor entry (binding)", loc="left", fontsize=9.5)
    fig.suptitle("Choice of the transfer-time split (green: both limits met; the cheapest such split is chosen)",
                 x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def cam_figure(n, plan, cam_rows, path: Path, dv_cam, cw_mod):
    S.apply()
    fig, ax = plt.subplots(figsize=(10, 4.6))
    sel = np.where(np.isin(plan.phase_of_t, ["leg 2 corridor", "hold (go/no-go)", "final translation"]))[0][::25]
    T = 2 * np.pi / n * 2
    tt = np.arange(0, T, 5.0)
    for i in sel:
        X = cw_mod.propagate(n, np.r_[plan.r[i], plan.v[i] + dv_cam * np.array([0, 1.0, 0])], tt)
        ax.plot(X[:, 1], X[:, 0], color=S.GREEN, lw=0.6, alpha=0.7)
    for ph, col in PHASE_COLORS.items():
        k = plan.phase_of_t == ph
        if k.any():
            ax.plot(plan.r[k, 1], plan.r[k, 0], color=col, lw=2.0)
    ax.plot(0, 0, "s", color=S.INK)
    ax.set_xlabel("y along-track [m]")
    ax.set_ylabel("x radial [m]")
    ok = [r for r in cam_rows if r["cam_direction"].startswith("+Y")][0]
    ax.set_title(f"Breakout/CAM: +{dv_cam * 100:.0f} cm/s posigrade impulse from any point of corridor, hold or final "
                 f"translation\nnever closes further and stays >= {ok['closest_return_10min_to_2orbits_m']:.0f} m away after 10 min "
                 f"(2-orbit check)", loc="left", fontsize=9)
    ax.set_aspect("equal")
    fig.savefig(path)
    plt.close(fig)
    return path
