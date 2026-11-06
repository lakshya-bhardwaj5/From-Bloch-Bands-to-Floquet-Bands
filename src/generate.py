"""Generate Chapter 2 figures, data, and native GitHub README GIFs.

Run from any directory: python path/to/src/generate.py
No video encoder, display server, notebook, or network connection is needed.
"""

import argparse
import csv
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from PIL import Image
import scipy

from lattice import band_metrics, dispersion, eigensystem, fourier_hopping

ROOT = Path(__file__).resolve().parents[1]
ASSETS, DATA = ROOT / "assets", ROOT / "data"
TEAL, INK, GRAY, RULE = "#1F6F5F", "#171717", "#717773", "#E4E7E5"
COLORS = [TEAL, INK, GRAY]
STYLES = ["-", "--", ":"]
KAPPA = np.linspace(-1, 1, 401)
DEPTHS = [2, 4, 8, 12]


def theme():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 12,
        "axes.labelsize": 13, "axes.titlesize": 14,
        "xtick.labelsize": 11, "ytick.labelsize": 11,
        "text.color": INK, "axes.labelcolor": INK,
        "axes.edgecolor": GRAY, "axes.spines.top": False,
        "axes.spines.right": False, "axes.linewidth": .75,
        "legend.frameon": False, "legend.fontsize": 11,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "svg.fonttype": "none", "savefig.facecolor": "white",
    })


def k_axis(ax):
    ax.set_xlim(-1, 1)
    ax.set_xticks([-1, -.5, 0, .5, 1])
    ax.set_xlabel(r"Quasimomentum $k/k_L$")
    ax.set_ylabel(r"Energy $E_n/E_R$")
    ax.yaxis.set_major_locator(MaxNLocator(5))


def band_lines(ax, values, labels=False):
    return [ax.plot(KAPPA, values[:, n], color=COLORS[n], ls=STYLES[n],
                    lw=2.4 if n == 0 else 1.8,
                    label=f"Band n = {n}" if labels else None)[0]
            for n in range(3)]


def save_static(fig, stem):
    for ext in ["svg", "png"]:
        fig.savefig(ASSETS / f"{stem}.{ext}", dpi=180)
    plt.close(fig)


def static_figures():
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.0), sharex=True, sharey=True)
    fig.subplots_adjust(left=.075, right=.97, bottom=.09, top=.82,
                        hspace=.38, wspace=.23)
    fig.suptitle("Deeper lattices narrow bands and open gaps", x=.075,
                 ha="left", y=.975, fontsize=19, weight="semibold")
    fig.text(.075, .922, r"Static continuum lattice $V(x)=sE_R\sin^2(k_Lx)$; 25 plane waves",
             color=GRAY)
    arrays = []
    for ax, s in zip(axes.flat, DEPTHS):
        values = dispersion(KAPPA, s)
        arrays.append(values)
        band_lines(ax, values, labels=s == 2)
        met = band_metrics(s)
        ax.axhspan(values[-1, 0], values[-1, 1], color=TEAL, alpha=.075, lw=0)
        ax.set_title(f"s = {s}     lowest bandwidth = {met['bandwidth_ER']:.3f} E_R",
                     loc="left", fontsize=12, pad=9)
        k_axis(ax)
        ax.set_ylim(0, 16)
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="upper left",
               bbox_to_anchor=(.065, .894), ncol=3)
    np.savez_compressed(DATA / "bands.npz", kappa=KAPPA,
                        depths=DEPTHS, energies_ER=np.array(arrays))
    save_static(fig, "01_bands")

    depths = np.linspace(0, 12, 121)
    metrics = [band_metrics(float(s)) for s in depths]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0))
    fig.subplots_adjust(left=.08, right=.97, bottom=.18, top=.77, wspace=.3)
    fig.suptitle("The lowest band narrows as the first gap grows", x=.08,
                 ha="left", y=.97, fontsize=18, weight="semibold")
    fig.text(.08, .89, "Same absolute energy reference as Chapter 2; no fitted parameters", color=GRAY)
    axes[0].plot(depths, [r["bandwidth_ER"] for r in metrics], color=TEAL, lw=2.4)
    axes[0].set_ylabel(r"Lowest-band width $W_0/E_R$")
    axes[0].set_ylim(0, 1.05)
    axes[1].plot(depths, [r["edge_gap_ER"] for r in metrics], color=TEAL, lw=2.4,
                 label="25-plane-wave calculation")
    axes[1].plot(depths, depths / 2, color=GRAY, ls="--", lw=1.7,
                 label=r"Weak-lattice limit $s/2$")
    axes[1].set_ylabel(r"First edge gap $\Delta_{01}^{\rm edge}/E_R$")
    axes[1].set_ylim(0, 6.4)
    axes[1].legend(loc="upper left", fontsize=10)
    for ax in axes:
        ax.set_xlabel(r"Lattice depth $s=V_0/E_R$")
        ax.set_xlim(0, 12)
        ax.set_xticks([0, 4, 8, 12])
        ax.yaxis.set_major_locator(MaxNLocator(5))
    save_static(fig, "02_width_and_gap")
    with (DATA / "depth_sweep.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0), sharey=True)
    fig.subplots_adjust(left=.09, right=.97, bottom=.18, top=.74, wspace=.25)
    fig.suptitle("Nearest-neighbor bands improve in deeper lattices", x=.09,
                 ha="left", y=.97, fontsize=18, weight="semibold")
    fig.text(.09, .90, "Subtract the band mean and divide by the extracted hopping J1", color=GRAY)
    summary = []
    for ax, s in zip(axes, [2, 8]):
        energy = dispersion(KAPPA, s, bands=1)[:, 0]
        mean, hop = fourier_hopping(s)
        nearest = -2 * np.cos(np.pi * KAPPA)
        ax.plot(KAPPA, (energy - mean) / hop[0], color=TEAL, lw=2.4, label="Continuum band")
        ax.plot(KAPPA, nearest, color=INK, ls="--", lw=1.6, label="Nearest neighbor")
        err = np.max(np.abs(energy - (mean + hop[0] * nearest))) / np.ptp(energy)
        ax.set_title(f"s = {s}     max error / bandwidth = {100 * err:.2f}%",
                     loc="left", fontsize=12, pad=12)
        ax.set_xlim(-1, 1)
        ax.set_xticks([-1, -.5, 0, .5, 1])
        ax.set_xlabel(r"Quasimomentum $k/k_L$")
        ax.set_ylim(-2.4, 2.55)
        row = band_metrics(s)
        row.update(mean_ER=mean, J1_ER=float(hop[0]),
                   max_NN_error_over_bandwidth=float(err))
        summary.append(row)
    axes[0].set_ylabel(r"Scaled energy $(E_0-\overline{E})/J_1$")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper left",
               bbox_to_anchor=(.08, .858), ncol=2)
    save_static(fig, "03_tight_binding")
    # Full precision for the report's four reference depths.
    table = []
    for s in DEPTHS:
        row = band_metrics(s)
        row["J1_ER"] = float(fourier_hopping(s)[1][0])
        table.append(row)
    with (DATA / "chapter2_table.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    (DATA / "tight_binding.json").write_text(json.dumps(summary, indent=2) + "\n")


def capture(fig):
    fig.canvas.draw()
    return Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())


def save_gif(frames, stem, duration):
    """One shared palette keeps white backgrounds and line colors stable."""
    reference = frames[len(frames) // 2].quantize(colors=128)
    indexed = [frame.quantize(palette=reference, dither=Image.Dither.NONE)
               for frame in frames]
    indexed[0].save(ASSETS / f"{stem}.gif", save_all=True,
                    append_images=indexed[1:], duration=duration,
                    loop=0, optimize=True, disposal=1)


def depth_animation():
    fig, (potential_ax, bands_ax) = plt.subplots(
        1, 2, figsize=(11.6, 5.7), dpi=100, gridspec_kw={"width_ratios": [1, 1.35]})
    fig.subplots_adjust(left=.075, right=.975, bottom=.17, top=.74, wspace=.27)
    fig.text(.075, .955, "From free-particle crossings to Bloch bands", fontsize=20,
             weight="semibold")
    status = fig.text(.075, .875, "", fontsize=13, color=TEAL)
    fig.text(.075, .035, "Parameter sweep of stationary eigenproblems; animation frames are not time evolution.",
             fontsize=11, color=GRAY)
    x = np.linspace(-1.5, 1.5, 501)
    potential, = potential_ax.plot(x, np.zeros_like(x), color=INK, lw=2)
    potential_ax.set(xlim=(-1.5, 1.5), ylim=(0, 16), xlabel=r"Position $x/d$",
                     ylabel=r"Potential $V(x)/E_R$")
    potential_ax.set_xticks([-1, 0, 1])
    potential_ax.set_title("Standing-wave lattice", loc="left", fontsize=13)
    lines = band_lines(bands_ax, dispersion(KAPPA, 0), labels=True)
    k_axis(bands_ax)
    bands_ax.set_ylim(0, 16)
    bands_ax.set_title("Three lowest energy bands", loc="left", fontsize=13)
    fig.legend(*bands_ax.get_legend_handles_labels(), loc="upper left",
               bbox_to_anchor=(.064, .835), ncol=3, fontsize=11)
    # Smooth ping-pong loop, with pauses at the two endpoints.
    sequence = np.linspace(0, 12, 49)
    depths = np.r_[np.repeat(0., 6), sequence[1:], np.repeat(12., 6), sequence[-2:0:-1]]
    frames = []
    for s in depths:
        potential.set_ydata(s * np.sin(np.pi * x) ** 2)
        values = dispersion(KAPPA, float(s))
        for n, line in enumerate(lines):
            line.set_ydata(values[:, n])
        met = band_metrics(float(s))
        status.set_text(f"Depth s = {s:5.2f}     |     W0/E_R = {met['bandwidth_ER']:.3f}"
                        f"     |     edge gap / E_R = {met['edge_gap_ER']:.3f}")
        frames.append(capture(fig))
    save_gif(frames, "band_formation", 85)
    plt.close(fig)


def momentum_animation():
    s = 4.0
    values = dispersion(KAPPA, s)
    fig, (bands_ax, weights_ax) = plt.subplots(1, 2, figsize=(11.6, 5.7), dpi=100)
    fig.subplots_adjust(left=.075, right=.975, bottom=.17, top=.74, wspace=.27)
    fig.text(.075, .955, "Each point on a band is a different Bloch eigenstate",
             fontsize=19, weight="semibold")
    status = fig.text(.075, .875, "", fontsize=13, color=TEAL)
    fig.text(.075, .035, "Static lattice s = 4; reciprocal basis exp[i(k - 2 ell k_L)x]. No external force is applied.",
             fontsize=10.5, color=GRAY)
    band_lines(bands_ax, values, labels=True)
    k_axis(bands_ax)
    bands_ax.set_ylim(0, 12)
    guide = bands_ax.axvline(-1, color=GRAY, lw=.8, ls="--")
    dots = [bands_ax.plot([], [], "o", color=COLORS[n], ms=7)[0] for n in range(3)]
    fig.legend(*bands_ax.get_legend_handles_labels(), loc="upper left",
               bbox_to_anchor=(.064, .835), ncol=3, fontsize=11)
    ell = np.arange(-12, 13)
    show = np.abs(ell) <= 4
    bars = weights_ax.bar(ell[show], np.zeros(show.sum()), color=TEAL, width=.65)
    weights_ax.set(xlim=(-4.7, 4.7), ylim=(0, 1.02),
                   xlabel=r"Reciprocal-vector index $\ell$",
                   ylabel=r"Lowest-band weight $|a_{0\ell}|^2$")
    weights_ax.set_xticks([-4, -2, 0, 2, 4])
    weights_ax.set_yticks([0, .5, 1])
    weights_ax.set_title("Plane-wave components of band n = 0", loc="left", fontsize=12)
    sequence = np.linspace(-1, 1, 61)
    momenta = np.r_[sequence, sequence[-2:0:-1]]
    frames = []
    for k in momenta:
        energies, vectors = eigensystem(float(k), s, vectors=True)
        weights = np.abs(vectors[:, 0]) ** 2
        guide.set_xdata([k, k])
        for n, dot in enumerate(dots):
            dot.set_data([k], [energies[n]])
        for bar, weight in zip(bars, weights[show]):
            bar.set_height(weight)
        status.set_text(f"k/k_L = {k:+.2f}     |     E0/E_R = {energies[0]:.3f}"
                        f"     |     weight outside displayed ell: {weights[~show].sum():.1e}")
        frames.append(capture(fig))
    save_gif(frames, "quasimomentum_scan", 70)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--static-only", action="store_true", help="Skip both GIFs.")
    args = parser.parse_args()
    ASSETS.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    theme()
    static_figures()
    if not args.static_only:
        depth_animation()
        momentum_animation()
    (DATA / "environment.json").write_text(json.dumps({
        "python": platform.python_version(), "numpy": np.__version__,
        "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
        "pillow": Image.__version__, "cutoff_M": 12,
        "plane_waves": 25, "plot_momentum_points": 401,
        "fourier_momentum_points": 2048, "randomness": False,
    }, indent=2) + "\n")
    print("Generated static figures and data" + ("." if args.static_only else ", plus two looping GIFs."))


if __name__ == "__main__":
    main()
