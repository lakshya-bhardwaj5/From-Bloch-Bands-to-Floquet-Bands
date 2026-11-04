"""Independent limits, discretization refinements, and finite-domain checks.

This script raises AssertionError on a failed numerical acceptance condition.
Tolerances concern these specified finite calculations, not universal accuracy.
"""
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import mathieu_a, mathieu_b, j0, jv
from common import setup, static, band_fourier, dispersion, write_json
from single_band import KSTAR, numerical_average
from wavepackets import propagate, moments, integrated_hopping, PERIOD
from multiband import monodromy, strobe_probabilities

def main():
    setup();checks={}
    qgrid=np.linspace(-np.pi,np.pi,33)
    # Lowest THREE bands at all listed depths and momenta: M=6,8,10 vs M=12.
    static_refinement=[]
    for m in [6,8,10]:
        error=max(np.max(np.abs(static(q,s,m)[:3]-static(q,s,12)[:3])) for s in [2,4,8,12] for q in qgrid)
        static_refinement.append({"cutoff":m,"max_low3_error_ER":float(error)})
    assert static_refinement[-1]["max_low3_error_ER"]<1e-10
    checks["static_cutoff"]=static_refinement
    # Independent special-function benchmark at exact Bloch band edges.
    err=0.
    for s in [2,4,8,12]:
        expected=[s/2+mathieu_a(0,s/4),s/2+mathieu_b(1,s/4)]
        computed=[static(0,s)[0],static(np.pi,s)[0]]
        err=max(err,float(np.max(np.abs(np.array(expected)-computed))))
    assert err<1e-11;checks["mathieu_edge_error_ER"]=err
    # Free particle: folded plane waves, including degeneracies at q=0/pi.
    err=max(np.max(np.abs(static(q,0,8)[:5]-np.sort((q/np.pi+2*np.arange(-8,9))**2)[:5])) for q in qgrid)
    assert err<1e-12;checks["free_particle_error_ER"]=float(err)
    # Fourier representation: different momentum and hopping cutoffs.
    fourier=[]
    for s in [2,8]:
        q=np.linspace(-np.pi,np.pi,401)
        exact=np.array([static(x,s)[0] for x in q])
        for nq,nhop in [(512,10),(1024,20),(2048,20)]:
            _,_,mean,hops=band_fourier(s,nq=nq,nhop=nhop)
            err=np.max(np.abs(dispersion(q,mean,hops)-exact))
            fourier.append({"s":s,"nq":nq,"nhop":nhop,"max_error_ER":float(err)})
    assert max(r["max_error_ER"] for r in fourier if r["nhop"]==20)<1e-6
    checks["fourier_reconstruction"]=fourier
    # Period quadrature vs Bessel, and B(t) vs adaptive quadrature.
    err=max(abs(numerical_average(q,k,512)+2*j0(k)*np.cos(q)) for q in [-2.3,0,.8] for k in [0,1.5,KSTAR,3.5,8])
    assert err<1e-12;checks["period_average_error_J"]=float(err)
    err=0.
    for k in [0,1.5,KSTAR,3.5]:
        for t in [.173,.91,PERIOD,3.1*PERIOD]:
            direct=quad(lambda u:np.cos(k*np.sin(4*u)),0,t,epsabs=1e-12)[0]+1j*quad(lambda u:np.sin(k*np.sin(4*u)),0,t,epsabs=1e-12)[0]
            err=max(err,abs(integrated_hopping([t],k)[0]-direct))
    assert err<1e-11;checks["integrated_hopping_error_tauJ"]=float(err)
    harmonic_error=max(np.max(np.abs(integrated_hopping(np.linspace(0,20*PERIOD,161),k,nfourier=40)-integrated_hopping(np.linspace(0,20*PERIOD,161),k,nfourier=60))) for k in [0,1.5,KSTAR,3.5])
    assert harmonic_error<1e-12
    checks["integrated_hopping_40_vs_60_harmonics"]=float(harmonic_error)
    # Free NN spreading: P_j(t)=J_j(2*tau_J)^2, including signed j.
    nsite=512;j=np.arange(nsite)-nsite//2
    initial=np.zeros(nsite,complex);initial[nsite//2]=1
    t=np.array([0.,.23,2.7,13.])
    computed=np.abs(propagate(initial,t,0))**2
    exact=jv(j[None,:],2*t[:,None])**2
    err=np.max(np.abs(computed-exact));assert err<1e-12
    checks["Bessel_probability_error"]=float(err)
    # Orthogonal implementation: integrate nearest-neighbor real-space ODE.
    # On a ring, (Hc)_j=-exp(i*chi)c_(j+1)-exp(-i*chi)c_(j-1).
    def rhs(t,c):
        phase=np.exp(1j*1.5*np.sin(4*t))
        return 1j*(phase*np.roll(c,-1)+phase.conjugate()*np.roll(c,1))
    sol=solve_ivp(rhs,(0,3*PERIOD),initial,method="DOP853",rtol=1e-10,atol=1e-12,max_step=PERIOD/80)
    fft_final=propagate(initial,[3*PERIOD],1.5)[0]
    err=np.max(np.abs(sol.y[:,-1]-fft_final));assert sol.success and err<1e-9
    checks["real_space_ODE_vs_FFT_error"]=float(err)
    # Check width, norm, revival and ring size over the full plotted interval.
    finite=[]
    for n in [256,512,1024]:
        c=np.zeros(n,complex);c[n//2]=1
        times=np.linspace(0,20*PERIOD,161)
        states=propagate(c,times,0.)
        center,width,norm=moments(states)
        err=np.max(np.abs(width-np.sqrt(2)*times))
        finite.append({"Nsite":n,"max_width_error_sites":float(err),"max_norm_error":float(np.max(abs(norm-1)))})
    assert finite[-1]["max_width_error_sites"]<1e-9
    checks["ring_size"]=finite
    c=np.zeros(512,complex);c[256]=1
    revived=propagate(c,np.arange(21)*PERIOD,KSTAR)
    checks["collapse_revival_max_amplitude_error"]=float(np.max(np.abs(revived-c)))
    assert checks["collapse_revival_max_amplitude_error"]<1e-12
    # Continuum refinements focus on the most populated resonant case, plus
    # the long off-resonant period; error refers to the entire m=0,...,20 trace.
    s=8.;k=.35
    e,v=static(0,s,8,vectors=True);gap=e[1]-e[0]
    refinement=[]
    for om in [1.,gap]:
        reference=strobe_probabilities(monodromy(0,s,k,om,8,4096),v[:,0])
        for steps in [512,1024,2048]:
            u=monodromy(0,s,k,om,8,steps)
            prob=strobe_probabilities(u,v[:,0])
            refinement.append({"Omega_R":float(om),"nstep":steps,"max_Pexc_error_vs_4096":float(np.max(np.abs(prob-reference))),"unitarity_error":float(np.linalg.norm(u.conj().T@u-np.eye(len(u)),2))})
    assert max(r["max_Pexc_error_vs_4096"] for r in refinement if r["nstep"]==2048)<1e-4
    checks["continuum_time_refinement"]=refinement
    cutoff_check=[]
    for om in [1.,gap]:
        traces=[]
        for m in [6,8,10]:
            en,vec=static(0,s,m,vectors=True)
            traces.append(strobe_probabilities(monodromy(0,s,k,om,m,2048),vec[:,0]))
        cutoff_check.append({"Omega_R":float(om),"M6_vs_M10":float(np.max(np.abs(traces[0]-traces[2]))),"M8_vs_M10":float(np.max(np.abs(traces[1]-traces[2])))})
    assert max(r["M8_vs_M10"] for r in cutoff_check)<1e-7
    checks["continuum_plane_wave_refinement"]=cutoff_check
    # Check zero shaking monodromy against exact time-independent evolution.
    u=monodromy(0,s,0.,1.,8,128)
    exact=(v*np.exp(-2j*np.pi*e))@v.T
    err=np.linalg.norm(u-exact,2);assert err<1e-10
    checks["undriven_monodromy_error"]=float(err)
    checks["status"]="passed"
    checks["scope"]="Only the specified depths, momenta, cutoffs, times and drive strengths; no proof of convergence outside them."
    write_json("convergence.json",checks)
    print("All numerical acceptance checks passed.")

if __name__=="__main__":main()
