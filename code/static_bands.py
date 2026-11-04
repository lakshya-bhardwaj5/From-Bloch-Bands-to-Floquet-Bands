"""Figures 1-2, static-band data and a numerical tight-binding assessment."""
import numpy as np
import matplotlib.pyplot as plt
from common import setup, static, band_fourier, dispersion, savefig, write_json, DATA, TEAL, GRAY, INK

def main():
    setup()
    q = np.linspace(-np.pi, np.pi, 401)
    fig, axes = plt.subplots(2, 2, figsize=(6.6, 4.9), sharex=True, layout="constrained")
    rows = []
    archive = {"q": q}
    for ax, s in zip(axes.flat, [2, 4, 8, 12]):
        e = np.array([static(x, s)[:3] for x in q])
        archive[f"s{s}"] = e
        for n, color in enumerate([TEAL, GRAY, INK]):
            ax.plot(q/np.pi, e[:, n], color=color, ls=["-", "--", ":"][n], label=f"n={n}")
        ax.set(title=f"Depth s = {s}", ylabel=r"Energy $E_n/E_R$", xticks=[-1, 0, 1])
        _, en, mean, hops = band_fourier(s)
        approx = mean - 2*hops[0]*np.cos(q)
        nn_error = np.max(np.abs(e[:, 0] - approx))
        rows.append({"s": s, "E0_center": float(e[200, 0]), "E0_edge": float(e[0, 0]),
                     "W0": float(np.ptp(en)), "gap_edge": float(e[0, 1]-e[0, 0]),
                     "gap_center": float(e[200, 1]-e[200, 0]), "Ebar": mean,
                     "J1": float(hops[0]), "J2": float(hops[1]),
                     "J3": float(hops[2]), "nn_max_error": float(nn_error),
                     "nn_relative_error": float(nn_error/np.ptp(en))})
    axes[0, 0].legend(loc="upper center", ncol=3)
    for ax in axes[-1]: ax.set_xlabel(r"Quasimomentum $k/k_L=q/\pi$")
    savefig(fig, "01_static_bands")
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.75), layout="constrained")
    for ax, s in zip(axes, [2, 8]):
        qs, energy, mean, hops = band_fourier(s)
        ax.plot(qs/np.pi, (energy-mean)/hops[0], color=TEAL, label="Plane-wave band")
        ax.plot(qs/np.pi, -2*np.cos(qs), color=INK, ls="--", label="Nearest neighbor")
        ax.set(title=f"Depth s = {s}", xlabel=r"$k/k_L$", ylabel=r"$(E_0-\bar E)/J_1$", xticks=[-1, 0, 1])
    axes[0].legend(loc="upper center")
    savefig(fig, "02_tight_binding")
    np.savez_compressed(DATA / "static_bands.npz", **archive)
    write_json("static_summary.json", rows)
    np.savetxt(DATA / "static_summary.csv", [[r[x] for x in ["s","E0_center","E0_edge","W0","gap_edge","J1","J2","nn_relative_error"]] for r in rows],
               delimiter=",", header="s,E0_center_ER,E0_edge_ER,W0_ER,gap_edge_ER,J1_ER,J2_ER,max_NN_error_over_W0", comments="")

if __name__ == "__main__": main()
