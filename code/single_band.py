"""Exact nearest-neighbor Floquet bands and isolated-band longer-range correction."""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0, jn_zeros
from common import setup, band_fourier, dispersion, savefig, write_json, DATA, TEAL, INK, GRAY

KSTAR = float(jn_zeros(0, 1)[0])
KVALUES = [0., 1.5, KSTAR, 3.5]

def numerical_average(q, strength, ntime=2048):
    """Periodic trapezoid quadrature, excluding the duplicated endpoint; units J."""
    phase = 2*np.pi*np.arange(ntime)/ntime
    return np.mean(-2*np.cos(np.asarray(q)[..., None] + strength*np.sin(phase)), axis=-1)

def main():
    setup()
    q = np.linspace(-np.pi, np.pi, 401)
    fig, axes = plt.subplots(2, 2, figsize=(6.6, 4.5), sharex=True, sharey=True, layout="constrained")
    errors = []
    for ax, strength in zip(axes.flat, KVALUES):
        ax.axhline(0, color="#D8D8D4", lw=.6)
        ax.plot(q/np.pi, -2*j0(strength)*np.cos(q), color=TEAL, label="Bessel formula")
        points = q[::20]
        avg = numerical_average(points, strength)
        ax.plot(points/np.pi, avg, "o", color=INK, mfc="white", ms=3, mew=.7, label="Period quadrature")
        errors.append(float(np.max(np.abs(avg+2*j0(strength)*np.cos(points)))))
        ax.set(title=f"K = {strength:.5g}", ylabel=r"Quasienergy $(\varepsilon-\bar E)/J$", xticks=[-1, 0, 1], ylim=(-2.2,2.2))
    axes[0,0].legend(loc="upper center", fontsize=8)
    for ax in axes[-1]: ax.set_xlabel(r"Quasimomentum $k/k_L$")
    savefig(fig, "03_floquet_bands")
    # Include the roots themselves so a connected line does not miss a cusp.
    kgrid = np.unique(np.r_[np.linspace(0, 8, 1001), jn_zeros(0, 2)])
    fig, axes = plt.subplots(2,1,figsize=(6.6,4.4),sharex=True,layout="constrained")
    axes[0].plot(kgrid, j0(kgrid), color=TEAL)
    axes[0].axhline(0,color=INK,lw=.6)
    axes[0].set(ylabel=r"Signed tunneling $J_{\rm eff}/J$")
    axes[1].plot(kgrid, np.abs(j0(kgrid)), color=TEAL)
    axes[1].set(xlabel=r"Driving strength $K=F_0d/(\hbar\omega)$", ylabel=r"Bandwidth $W/(4J)$")
    for ax in axes:
        for root in jn_zeros(0,2): ax.axvline(root,color=GRAY,ls=":",lw=.8)
        ax.set_xlim(0,8)
    savefig(fig,"04_tunneling_bandwidth")
    # Retain 20 Fourier hoppings: this is an isolated-band calculation, not a
    # full continuum Floquet calculation. Interband transitions are omitted.
    fig, axes = plt.subplots(1,2,figsize=(6.6,2.85),layout="constrained")
    collapse = []
    archive = {"K":kgrid}
    qdense = np.linspace(-np.pi,np.pi,2048,endpoint=False)
    for ax,s in zip(axes,[2,8]):
        _,energy,mean,hops = band_fourier(s)
        l=np.arange(1,len(hops)+1)
        basis=np.cos(qdense[:,None]*l)
        widths=np.array([np.ptp(-2*basis@(hops*j0(l*k))) for k in kgrid])
        w0=np.ptp(energy)
        ax.plot(kgrid,widths/w0,color=TEAL,label="20 hopping ranges")
        ax.plot(kgrid,np.abs(j0(kgrid)),color=INK,ls="--",label="Nearest neighbor")
        ax.set(title=f"Depth s = {s}",xlabel=r"Driving strength $K$",ylabel=r"Relative bandwidth $W(K)/W(0)$",xlim=(1.9,2.9),ylim=(0,.3))
        residual=np.ptp(-2*basis@(hops*j0(l*KSTAR)))
        ax.axvline(KSTAR,color=GRAY,ls=":",lw=.8)
        collapse.append({"s":s,"W_at_Kstar_ER":float(residual),"W_at_Kstar_over_W0":float(residual/w0),"J2_ER":float(hops[1])})
        archive[f"width_s{s}"]=widths
        if s==8:
            zoomk=np.linspace(KSTAR-.02,KSTAR+.02,241)
            zoomw=np.array([np.ptp(-2*basis@(hops*j0(l*k))) for k in zoomk])/w0
            inset=ax.inset_axes([.34,.57,.42,.37])
            inset.plot(zoomk,zoomw,color=TEAL,lw=1.1)
            inset.plot(zoomk,np.abs(j0(zoomk)),color=INK,ls="--",lw=.9)
            inset.axvline(KSTAR,color=GRAY,ls=":",lw=.6)
            inset.set(xlim=(KSTAR-.02,KSTAR+.02),ylim=(0,.013),xticks=[2.39,2.42],yticks=[0,.01])
            inset.tick_params(labelsize=8)
    axes[0].legend(loc="upper left",fontsize=8)
    savefig(fig,"09_long_range_collapse")
    np.savez_compressed(DATA/"single_band.npz",**archive)
    write_json("single_band_summary.json",{"Kstar":KSTAR,"strengths":KVALUES,"J0_values":[float(j0(k)) for k in KVALUES],"quadrature_max_errors":errors,"long_range":collapse})

if __name__=="__main__": main()
