"""
3-Axis Integrative Framework (3-AIF) ODE Simulation
====================================================
A six-variable ODE system modelling difficult-to-treat (D2T) rheumatic disease.

State variables (normalised to [0, 1]):
    T - Tolerance (Axis 1), D - Danger signal (Axis 2), S - Stress (Axis 3),
    E - Energy reserve, R - Recovery capacity, M - Microbiota diversity

Reference: Jung S. bioRxiv (2026).
Author: Sungsoo Jung, MD PhD — SoonChunHyang University Bucheon Hospital
ORCID: 0009-0006-6885-192X
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import os

# === Default parameters (Supplementary Table S4) ===
DEFAULT_PARAMS = {
    'alpha1': 0.55, 'beta1': 0.22, 'gamma1': 0.08, 'epsilon': 0.55,
    'alpha2': 0.18, 'beta2': 0.50, 'beta2env': 0.35, 'gamma2': 0.30,
    'alpha3': 0.40, 'beta3': 0.55, 'gamma3': 0.04,
    'delta1': 0.65, 'delta2': 0.30, 'delta3': 0.15, 'delta4': 0.15,
    'epsilon1': 0.40, 'epsilon2': 0.30, 'epsilon3': 0.25,
    'mu1': 0.45, 'mu2': 0.25, 'mu3': 0.08,
    'Km': 0.50, 'Ka': 0.50,
}
INMAS_MOD = {'alpha2_f': 0.45, 'delta1_f': 0.50, 'beta3_f': 2.0, 'env_f': 0.20}
HEALTHY_IC = [0.90, 0.05, 0.05, 0.80, 0.75, 0.85]

def mTORC1(E, Km=0.5): return E**2 / (Km**2 + E**2)
def AMPK(E, Ka=0.5):   return (1-E)**2 / (Ka**2 + (1-E)**2)

def ode_system(t, y, params, env_func, intv_t=None):
    T,D,S,E,R,M = [np.clip(v, 1e-10, 1-1e-10) for v in y]
    env = env_func(t)
    p = params; a2,d1,b3 = p['alpha2'],p['delta1'],p['beta3']
    if intv_t and t >= intv_t:
        a2 *= INMAS_MOD['alpha2_f']; d1 *= INMAS_MOD['delta1_f']
        b3 *= INMAS_MOD['beta3_f'];  env *= INMAS_MOD['env_f']
    mT, aK = mTORC1(E, p['Km']), AMPK(E, p['Ka'])
    return [
        p['alpha1']*M*(1-T) - p['beta1']*D*T - p['gamma1']*S*T - p['epsilon']*env*T,
        a2*(1-T)*E + p['beta2env']*env - p['beta2']*D + p['gamma2']*S*(1-R),
        p['alpha3']*D*(1-S) - b3*R*S + p['gamma3']*mT*(1-S),
        -d1*D*E + p['delta2']*R*(1-E) - p['delta3']*S*E + p['delta4']*(1-E),
        p['epsilon1']*aK*(1-R) - p['epsilon2']*S*R - p['epsilon3']*D*R,
        p['mu1']*T*(1-M) - p['mu2']*D*M - p['mu3']*(1-E)*M,
    ]

def env_healthy(t):     return 0.50 if 5<=t<=8 else 0.0
def env_chronic(t):     return 0.65 if t>=5 else 0.0
def env_therapeutic(t): return 0.65 if t>=5 else 0.0

def run_simulation(scenario='A', t_span=(0,80), t_points=500, params=None):
    if params is None: params = DEFAULT_PARAMS.copy()
    cfg = {'A':(env_healthy,None), 'B':(env_chronic,None), 'C':(env_therapeutic,30)}
    ef, it = cfg[scenario]
    return solve_ivp(lambda t,y: ode_system(t,y,params,ef,it), t_span, HEALTHY_IC,
                     t_eval=np.linspace(*t_span, t_points), method='RK45',
                     rtol=1e-8, atol=1e-10, max_step=0.1)

LABELS = ['T (Tolerance)','D (Danger)','S (Stress)','E (Energy)','R (Recovery)','M (Microbiota)']
COLORS = ['#2196F3','#F44336','#9C27B0','#FF9800','#4CAF50','#795548']

def plot_trajectories(save_path=None):
    """Figure 2: Three disease trajectory scenarios."""
    fig, axes = plt.subplots(1,3,figsize=(16,5))
    titles = ['A. Healthy Resolution','B. Chronic Disease (D2T)','C. Therapeutic Intervention']
    for i,sc in enumerate(['A','B','C']):
        sol = run_simulation(sc); ax = axes[i]
        for j,(lb,co) in enumerate(zip(LABELS,COLORS)):
            ax.plot(sol.t, sol.y[j], label=lb, color=co, linewidth=1.5)
        ax.set_xlabel('Time (a.u.)'); ax.set_ylabel('Value [0,1]')
        ax.set_title(titles[i], fontweight='bold'); ax.set_ylim(-0.05,1.05)
        ax.legend(fontsize=7, loc='right'); ax.grid(True, alpha=0.3)
        if sc=='C': ax.axvline(30, color='green', ls='--', alpha=0.6, lw=1.5)
    plt.tight_layout()
    if save_path: plt.savefig(save_path, dpi=300, bbox_inches='tight'); print(f"  Saved: {save_path}")
    plt.close()

def plot_metabolic_switch(save_path=None):
    """Figure 3B: mTORC1-AMPK crossover."""
    E = np.linspace(0,1,200)
    fig, ax = plt.subplots(figsize=(8,5))
    ax.plot(E, [mTORC1(e) for e in E], 'r-', lw=2.5, label='mTORC1(E)')
    ax.plot(E, [AMPK(e) for e in E], 'b-', lw=2.5, label='AMPK(E)')
    ax.axvline(0.5, color='gold', ls='--', lw=1.5, alpha=0.7, label='Crossover')
    ax.fill_between(E,0,1,where=E<0.5,alpha=0.08,color='blue')
    ax.fill_between(E,0,1,where=E>=0.5,alpha=0.08,color='red')
    ax.set_xlabel('Energy (E)'); ax.set_ylabel('Activity')
    ax.set_title('mTORC1-AMPK Metabolic Switch', fontweight='bold')
    ax.legend(); ax.grid(True, alpha=0.3); plt.tight_layout()
    if save_path: plt.savefig(save_path, dpi=300, bbox_inches='tight'); print(f"  Saved: {save_path}")
    plt.close()

def plot_bifurcation(save_path=None):
    """Figure 3C: Bifurcation diagram."""
    p = DEFAULT_PARAMS.copy(); envs = np.linspace(0,0.8,40)
    yf = list(HEALTHY_IC); Df, Db = [], []
    for ev in envs:
        s = solve_ivp(lambda t,y: ode_system(t,y,p,lambda t:ev),(0,200),yf,method='RK45',rtol=1e-8,atol=1e-10,max_step=0.5)
        yf = s.y[:,-1].tolist(); Df.append(yf[1])
    yb = list(yf)
    for ev in reversed(envs):
        s = solve_ivp(lambda t,y: ode_system(t,y,p,lambda t:ev),(0,200),yb,method='RK45',rtol=1e-8,atol=1e-10,max_step=0.5)
        yb = s.y[:,-1].tolist(); Db.append(yb[1])
    Db.reverse()
    fig, ax = plt.subplots(figsize=(8,5))
    ax.plot(envs,Df,'r-o',lw=2,ms=4,label='Forward',alpha=0.8)
    ax.plot(envs,Db,'b--s',lw=2,ms=4,label='Backward',alpha=0.8)
    ax.set_xlabel('Environmental Stress'); ax.set_ylabel('Steady-State D')
    ax.set_title('Bifurcation Diagram', fontweight='bold')
    ax.legend(); ax.grid(True,alpha=0.3); plt.tight_layout()
    if save_path: plt.savefig(save_path, dpi=300, bbox_inches='tight'); print(f"  Saved: {save_path}")
    plt.close()

def plot_sensitivity(variation=0.20, save_path=None):
    """Figure 4A: Sensitivity tornado plot."""
    pb = DEFAULT_PARAMS.copy()
    sb = run_simulation('B', params=pb)
    DI0 = 0.4*sb.y[1,-1] + 0.3*sb.y[2,-1] + 0.3*(1-sb.y[0,-1])
    test = ['alpha1','alpha2','alpha3','beta1','beta2','beta3','gamma1','gamma2','gamma3',
            'delta1','delta2','delta3','delta4','epsilon1','epsilon2','epsilon3','mu1','mu2','mu3']
    res = {}
    for pn in test:
        dv = []
        for f in [1+variation, 1-variation]:
            pm = pb.copy(); pm[pn] *= f; s = run_simulation('B', params=pm)
            dv.append(0.4*s.y[1,-1]+0.3*s.y[2,-1]+0.3*(1-s.y[0,-1]) - DI0)
        res[pn] = {'high':dv[0],'low':dv[1],'range':abs(dv[0]-dv[1])}
    sp = sorted(res, key=lambda x:res[x]['range'], reverse=True)
    fig, ax = plt.subplots(figsize=(10,8))
    for i,pn in enumerate(sp):
        r = res[pn]
        ax.barh(i, r['high'], height=0.4, color='#F44336' if r['high']>0 else '#4CAF50', alpha=0.8)
        ax.barh(i, r['low'],  height=0.4, color='#4CAF50' if r['low']<0 else '#F44336', alpha=0.8)
    ax.set_yticks(range(len(sp))); ax.set_yticklabels(sp)
    ax.set_xlabel('ΔDisease Index'); ax.set_title('Sensitivity: ±20% Variation', fontweight='bold')
    ax.axvline(0,color='k',lw=0.5); ax.grid(True,alpha=0.3,axis='x'); ax.invert_yaxis()
    for i in range(3): ax.get_yticklabels()[i].set_fontweight('bold'); ax.get_yticklabels()[i].set_color('#D32F2F')
    plt.tight_layout()
    if save_path: plt.savefig(save_path, dpi=300, bbox_inches='tight'); print(f"  Saved: {save_path}")
    plt.close()
    print("  Top 3:"); [print(f"    {p}: {res[p]['range']:.4f}") for p in sp[:3]]

if __name__ == '__main__':
    out = 'figures'; os.makedirs(out, exist_ok=True)
    print("="*60+"\n  3-AIF ODE Simulation\n"+"="*60)
    print("\n[0] Steady states:")
    for sc in 'ABC':
        v = run_simulation(sc).y[:,-1]
        print(f"  {sc}: T={v[0]:.3f} D={v[1]:.3f} S={v[2]:.3f} E={v[3]:.3f} R={v[4]:.3f} M={v[5]:.3f}")
    print("\n[1] Figure 2..."); plot_trajectories(f'{out}/fig2_trajectories.png')
    print("[2] Figure 3B..."); plot_metabolic_switch(f'{out}/fig3b_metabolic_switch.png')
    print("[3] Figure 3C..."); plot_bifurcation(f'{out}/fig3c_bifurcation.png')
    print("[4] Figure 4A..."); plot_sensitivity(save_path=f'{out}/fig4a_sensitivity.png')
    print("\n"+"="*60+f"\n  Done. Figures in ./{out}/\n"+"="*60)
