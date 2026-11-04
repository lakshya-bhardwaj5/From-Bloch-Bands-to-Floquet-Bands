"""Exact-model Fourier propagation, probability maps, widths and directed motion.

No effective Hamiltonian is used to generate the continuous-time curves.
The full time-dependent nearest-neighbor dispersion is integrated first.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, j0
from common import setup, savefig, write_json, DATA, TEAL, INK, GRAY, COLORS, STYLES, CMAP
from single_band import KSTAR, KVALUES

OMEGA_J=4.0  # hbar*omega/J, exact NN model needs no high-frequency expansion
PERIOD=2*np.pi/OMEGA_J

def integrated_hopping(tau, strength, omega=OMEGA_J, nfourier=40):
    """B(tau)=integral_0^tau exp[i*K*sin(omega*u)] du in tau=J*t/hbar.

    The Jacobi-Anger expansion gives a rapidly convergent complex series.
    expm1 prevents cancellation for small arguments.
    """
    tau=np.atleast_1d(tau)
    n=np.r_[np.arange(-nfourier,0),np.arange(1,nfourier+1)]
    oscillating=np.expm1(1j*tau[:,None]*omega*n)
    return j0(strength)*tau + oscillating @ (jv(n,strength)/(1j*n*omega))

def propagate(initial, tau, strength, omega=OMEGA_J):
    """Unitary FFT pair; j is centered, q uses NumPy FFT frequency ordering."""
    nsite=len(initial)
    q=2*np.pi*np.fft.fftfreq(nsite)
    amplitude=np.fft.fft(np.fft.ifftshift(initial),norm="ortho")
    b=integrated_hopping(tau,strength,omega)
    # Integral of H/J=-2*cos(q+K*sin) is -2*Re(exp(iq)*B).
    phase=np.exp(2j*np.real(b[:,None]*np.exp(1j*q)))
    return np.fft.fftshift(np.fft.ifft(amplitude[None,:]*phase,axis=1,norm="ortho"),axes=1)

def moments(states):
    j=np.arange(len(states[0]))-len(states[0])//2
    prob=np.abs(states)**2
    center=prob@j
    variance=np.maximum(0.,prob@(j*j)-center**2)
    return center,np.sqrt(variance),np.sum(prob,axis=1)

def main():
    setup()
    nsite=1024
    j=np.arange(nsite)-nsite//2
    initial=np.zeros(nsite,complex);initial[nsite//2]=1
    tau=np.linspace(0,20*PERIOD,1601)
    strobe=np.arange(21)*PERIOD
    fig,axes=plt.subplots(2,2,figsize=(6.6,4.8),sharex=True,sharey=True,layout="constrained")
    widths=[];errors=[];archive={"tau_J":tau,"site":j,"strobe_tau_J":strobe}
    for ax,k in zip(axes.flat,KVALUES):
        state=propagate(initial,tau,k)
        probability=np.abs(state)**2
        center,width,norm=moments(state)
        widths.append(width)
        mask=np.abs(j)<=70
        # Fixed probability scale across all panels preserves comparison.
        im=ax.pcolormesh(j[mask],tau/PERIOD,probability[:,mask],cmap=CMAP,vmin=0,vmax=1,shading="auto",rasterized=True)
        ax.set(title=f"K = {k:.5g}",ylabel=r"Time $t/T$",xlim=(-70,70),ylim=(0,20))
        b=integrated_hopping(tau,k)
        errors.append({"K":k,"max_norm_error":float(np.max(np.abs(norm-1))),
                       "max_center_error":float(np.max(np.abs(center))),
                       "max_width_error":float(np.max(np.abs(width-np.sqrt(2)*np.abs(b)))),
                       "boundary_probability":float(np.max(np.sum(probability[:,np.abs(j)>nsite//2-20],axis=1)))})
        archive[f"prob_K{k:.6f}"]=probability[:,mask]
    for ax in axes[-1]:ax.set_xlabel(r"Lattice site $j=x/d$")
    fig.colorbar(im,ax=axes.ravel().tolist(),label=r"Site probability $P_j(t)$",fraction=.035,pad=.02)
    savefig(fig,"05_probability_evolution")
    # A second, linear probability scale exposes the expanding low-density
    # wings. The identical data are used; the saturation threshold is explicit.
    fig,axes=plt.subplots(1,2,figsize=(6.6,2.8),layout="constrained")
    for ax,k in zip(axes,[0.,3.5]):
        p=np.abs(propagate(initial,tau,k))**2;mask=np.abs(j)<=70
        im=ax.pcolormesh(j[mask],tau/PERIOD,p[:,mask],cmap=CMAP,vmin=0,vmax=.05,shading="auto",rasterized=True)
        ax.set(title=f"K = {k:g}; color saturates at 0.05",xlabel=r"Site $j$",ylabel=r"Time $t/T$",xlim=(-70,70))
    fig.colorbar(im,ax=axes.tolist(),label=r"$P_j$ (clipped color scale)",fraction=.035,pad=.02)
    savefig(fig,"05b_expanding_wings")
    fig,axes=plt.subplots(1,2,figsize=(6.6,2.9),layout="constrained")
    for k,width,color,ls in zip(KVALUES,widths,COLORS,STYLES):
        axes[0].plot(tau/PERIOD,width,color=color,ls=ls,label=f"K={k:.5g}")
        prediction=np.sqrt(2)*abs(j0(k))*strobe
        axes[0].plot(strobe/PERIOD,prediction,"o",mfc="white",color=color,ms=2.8,mew=.7)
    axes[0].set(xlabel=r"Time $t/T$",ylabel=r"RMS width $\sigma_x/d$",xlim=(0,20))
    axes[0].legend(loc="upper left")
    mic_tau=np.linspace(0,2*PERIOD,801)
    micro=propagate(initial,mic_tau,KSTAR)
    _,micro_width,_=moments(micro)
    axes[1].plot(mic_tau/PERIOD,micro_width,color=TEAL)
    axes[1].plot([0,1,2],[0,0,0],"o",mfc="white",color=INK,ms=4)
    axes[1].set(title=r"At $K=2.4048256$",xlabel=r"Time $t/T$",ylabel=r"RMS width $\sigma_x/d$",xticks=[0,.5,1,1.5,2])
    savefig(fig,"06_width_micromotion")
    # Snapshots independently make the complete-cycle revival visible.
    sample=np.array([0,.25,.5,.75,1.])*PERIOD
    p=np.abs(propagate(initial,sample,KSTAR))**2
    fig,ax=plt.subplots(figsize=(6.6,2.75),layout="constrained")
    for idx,(phase,color,ls) in enumerate(zip([0,.25,.5,.75,1.],[INK,GRAY,TEAL,GRAY,INK],["-","--","-",":","--"])):
        ax.plot(j,p[idx],color=color,ls=ls,marker=["o","s","^","v","x"][idx],ms=4,label=f"t/T={phase:g}")
    ax.set(xlim=(-3,3),xticks=np.arange(-3,4),xlabel="Lattice site j",ylabel=r"Probability $P_j$",ylim=(-.02,1.05))
    ax.legend(ncol=5,loc="upper center",bbox_to_anchor=(.5,1.16))
    savefig(fig,"07_localization_snapshots")
    # Same initial momentum in every run: re-preparing a band minimum would
    # remove the intended sign-sensitive comparison.
    sigma=6.;q0=np.pi/3
    gaussian=np.exp(-j*j/(4*sigma*sigma)+1j*q0*j)
    gaussian/=np.linalg.norm(gaussian)
    aq=np.fft.fft(np.fft.ifftshift(gaussian),norm="ortho")
    q=2*np.pi*np.fft.fftfreq(nsite)
    sinmean=float(np.sum(np.abs(aq)**2*np.sin(q)))
    fig,axes=plt.subplots(1,2,figsize=(6.6,2.9),layout="constrained")
    drift=[]
    for k,color,ls in zip(KVALUES,COLORS,STYLES):
        center,width,norm=moments(propagate(gaussian,tau,k))
        cs,_,_=moments(propagate(gaussian,strobe,k))
        slope=2*j0(k)*sinmean
        axes[0].plot(tau/PERIOD,center,color=color,ls=ls,label=f"K={k:.5g}")
        axes[0].plot(strobe/PERIOD,slope*strobe,"o",mfc="white",color=color,ms=2.8,mew=.7)
        drift.append({"K":k,"predicted_drift_dJ_over_hbar":float(slope),
                      "fitted_drift_dJ_over_hbar":float(np.polyfit(strobe,cs,1)[0]),
                      "max_strobe_error_sites":float(np.max(np.abs(cs-slope*strobe)))})
    axes[0].set(xlabel=r"Time $t/T$",ylabel=r"Center $\langle x\rangle/d$",xlim=(0,20))
    axes[0].legend(loc="upper left")
    for k,color,ls in [(0.,INK,STYLES[0]),(3.5,COLORS[3],STYLES[3])]:
        p=np.abs(propagate(gaussian,[20*PERIOD],k)[0])**2
        axes[1].plot(j,p,color=color,ls=ls,label=f"K={k:g}")
    axes[1].set(xlabel=r"Site $j$",ylabel=r"$P_j(20T)$",xlim=(-45,75))
    axes[1].legend(loc="upper center")
    savefig(fig,"08_directed_motion")
    archive["widths"]=np.array(widths)
    np.savez_compressed(DATA/"wavepackets.npz",**archive)
    write_json("wavepacket_summary.json",{"Omega_J":OMEGA_J,"Nsite":nsite,"delta_checks":errors,"gaussian_sigma_sites":sigma,"q0":q0,"sin_q_mean":sinmean,"drift":drift,"max_collapse_micromotion_width_sites":float(np.max(micro_width))})

if __name__=="__main__":main()
