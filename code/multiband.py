"""Bounded continuum calculation of interband excitation at k=0.

Use the unitary exponential-midpoint method over ONE period. Powers of its
monodromy matrix then give observations at integer periods. This includes all
2*M+1 plane-wave states, with no projection onto the lowest static band.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal, schur
from common import setup, static, savefig, write_json, DATA, TEAL, INK, GRAY

def monodromy(q,s,strength,omega,cutoff=8,nstep=1024):
    """All energies in E_R; omega=hbar*angular_frequency/E_R; q=kd.

    U <- exp[-i H(t_mid)*dt] U. Every factor is unitary (to roundoff).
    The second-order time discretization is explicitly checked by refinement.
    """
    ell=np.arange(-cutoff,cutoff+1)
    p=q/np.pi+2*ell
    off=np.full(2*cutoff,-s/4)
    step=2*np.pi/(omega*nstep)
    u=np.eye(len(ell),dtype=complex)
    a=strength/np.pi
    for it in range(nstep):
        shift=a*np.sin(omega*(it+.5)*step)
        e,v=eigh_tridiagonal((p+shift)**2+s/2,off)
        u=v @ (np.exp(-1j*step*e)[:,None]*(v.T@u))
    return u

def strobe_probabilities(u,initial,nperiod=20):
    state=initial.astype(complex)
    p=[]
    for n in range(nperiod+1):
        p.append(1-abs(np.vdot(initial,state))**2)
        state=u@state
    return np.array(p)

def main():
    setup()
    s=8.;q=0.;strength=.35;cutoff=8;nstep=2048
    e,v=static(q,s,cutoff,vectors=True)
    gap=float(e[1]-e[0]);ground=v[:,0]
    frequencies=[1.,gap]
    pmat=2*np.arange(-cutoff,cutoff+1)
    p10=float(abs(np.vdot(v[:,1],pmat*ground)))
    fig,axes=plt.subplots(1,2,figsize=(6.6,3.1),layout="constrained")
    results=[];archive={"period":np.arange(21)}
    for ax,om in zip(axes,frequencies):
        u=monodromy(q,s,strength,om,cutoff,nstep)
        prob=strobe_probabilities(u,ground)
        ax.plot(np.arange(21),prob,"o-",color=TEAL,ms=3,label="Continuum")
        ax.axhline(0,color=INK,ls="--",lw=1,label="Isolated band")
        if om==gap:
            t=np.arange(21)*2*np.pi/om
            rwa=np.sin(strength/np.pi*p10*t)**2
            ax.plot(np.arange(21),rwa,color=GRAY,ls=":",label="Two-state RWA")
        # Schur decomposition of this nearly normal matrix gives stable,
        # orthonormal eigenvectors for eigenphase/overlap diagnostics.
        triangular,z=schur(u,output="complex")
        eps=-np.angle(np.diag(triangular))*om/(2*np.pi)
        weight=np.abs(z.conj().T@ground)**2
        results.append({"Omega_R":om,"max_Pexc":float(np.max(prob)),"Pexc_20T":float(prob[-1]),
                        "unitarity_error":float(np.linalg.norm(u.conj().T@u-np.eye(len(u)),2)),
                        "largest_static_ground_Floquet_weight":float(np.max(weight)),
                        "reduced_quasienergies_ER":eps.tolist(),"ground_weights":weight.tolist()})
        archive[f"Pexc_Omega{om:.6f}"]=prob
        ax.set(title=fr"$\hbar\omega/E_R={om:.5g}$",xlabel=r"Completed periods $m=t/T$",ylabel=r"Excited-band probability $P_{\rm ex}(mT)$",xlim=(0,20))
    handles,labels=axes[1].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside upper center",ncol=3,fontsize=8)
    savefig(fig,"10_multiband_excitation")
    np.savez_compressed(DATA/"multiband.npz",**archive)
    write_json("multiband_summary.json",{"s":s,"q":q,"K":strength,"cutoff":cutoff,"nstep":nstep,"gap_at_k0_ER":gap,"p10_hbar_kL":p10,"cases":results})

if __name__=="__main__":main()
