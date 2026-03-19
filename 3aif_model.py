import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# ============================================================
# Supplementary Table S4 default parameters
# ============================================================
DEFAULT_PARAMS = {
    "alpha1": 0.5,
    "alpha2": 0.8,
    "alpha3": 0.6,

    "beta1": 0.4,
    "beta2": 0.6,
    "beta2env": 0.5,
    "beta3": 0.5,

    "gamma1": 0.3,
    "gamma2": 0.4,
    "gamma3": 0.4,

    "delta1": 0.7,
    "delta2": 0.5,
    "delta3": 0.4,
    "delta4": 0.1,

    "eps": 0.2,
    "eps1": 0.6,
    "eps2": 0.4,
    "eps3": 0.3,

    "mu1": 0.5,
    "mu2": 0.4,
    "mu3": 0.3,

    "Km": 0.5,
    "Ka": 0.5,
}


# ============================================================
# Hill functions
# ============================================================
def mtorc1(E: float, Km: float) -> float:
    return E**2 / (Km**2 + E**2)


def ampk(E: float, Ka: float) -> float:
    return (1.0 - E)**2 / (Ka**2 + (1.0 - E)**2)


# ============================================================
# Environment functions by scenario
# Scenario A: env(t) = 0.5 for t in [5, 15], else 0
# Scenario B: env(t) = 0.7 for t in [0, 100]
# Scenario C: same as B
# ============================================================
def env_scenario_A(t: float) -> float:
    return 0.5 if 5.0 <= t <= 15.0 else 0.0


def env_scenario_B(t: float) -> float:
    return 0.7 if 0.0 <= t <= 100.0 else 0.0


def env_scenario_C(t: float) -> float:
    return 0.7 if 0.0 <= t <= 100.0 else 0.0


# ============================================================
# Time-dependent parameter logic for Scenario C
# Before t=30:
#   alpha2 = 0.8, delta1 = 0.7, beta3 = 0.5
# After t=30:
#   alpha2 = 0.6, delta1 = 0.4, beta3 = 1.2
# ============================================================
def scenario_C_params(t: float, base_params: dict) -> dict:
    p = dict(base_params)
    if t >= 30.0:
        p["alpha2"] = 0.6
        p["delta1"] = 0.4
        p["beta3"] = 1.2
    return p


# ============================================================
# ODE system
# State vector y = [T, D, S, E, R, M]
# ============================================================
def model_3aif(t: float, y: np.ndarray, params: dict, env_fn, param_schedule_fn=None):
    T, D, S, E, R, M = y

    # Optional numerical clipping because states are normalized to [0, 1]
    T = np.clip(T, 0.0, 1.0)
    D = np.clip(D, 0.0, 1.0)
    S = np.clip(S, 0.0, 1.0)
    E = np.clip(E, 0.0, 1.0)
    R = np.clip(R, 0.0, 1.0)
    M = np.clip(M, 0.0, 1.0)

    p = dict(params) if param_schedule_fn is None else param_schedule_fn(t, params)
    env_t = env_fn(t)

    mtor = mtorc1(E, p["Km"])
    amp = ampk(E, p["Ka"])

    dT = (
        p["alpha1"] * M * (1.0 - T)
        - p["beta1"] * D * T
        - p["gamma1"] * S * T
        - p["eps"] * env_t * T
    )

    dD = (
        p["alpha2"] * (1.0 - T) * E
        + p["beta2env"] * env_t
        - p["beta2"] * D
        + p["gamma2"] * S * (1.0 - R)
    )

    dS = (
        p["alpha3"] * D * (1.0 - S)
        - p["beta3"] * R * S
        + p["gamma3"] * mtor * (1.0 - S)
    )

    dE = (
        -p["delta1"] * D * E
        + p["delta2"] * R * (1.0 - E)
        - p["delta3"] * S * E
        + p["delta4"] * (1.0 - E)
    )

    dR = (
        p["eps1"] * amp * (1.0 - R)
        - p["eps2"] * S * R
        - p["eps3"] * D * R
    )

    dM = (
        p["mu1"] * T * (1.0 - M)
        - p["mu2"] * D * M
        - p["mu3"] * (1.0 - E) * M
    )

    return [dT, dD, dS, dE, dR, dM]


# ============================================================
# Simulation helper
# ============================================================
def run_simulation(
    scenario_name: str,
    y0: list,
    env_fn,
    params: dict = None,
    param_schedule_fn=None,
    t_span=(0.0, 100.0),
    n_eval=1000,
    method="RK45",
):
    if params is None:
        params = DEFAULT_PARAMS

    t_eval = np.linspace(t_span[0], t_span[1], n_eval)

    sol = solve_ivp(
        fun=lambda t, y: model_3aif(
            t=t,
            y=y,
            params=params,
            env_fn=env_fn,
            param_schedule_fn=param_schedule_fn,
        ),
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        method=method,
        vectorized=False,
    )

    if not sol.success:
        raise RuntimeError(f"{scenario_name} failed: {sol.message}")

    return sol


# ============================================================
# Scenario definitions from Supplementary Table S4
# ============================================================
SCENARIOS = {
    "A_Healthy_Resolution": {
        "y0": [0.8, 0.2, 0.2, 0.7, 0.6, 0.8],
        "env_fn": env_scenario_A,
        "param_schedule_fn": None,
    },
    "B_Chronic_D2T": {
        "y0": [0.5, 0.5, 0.4, 0.4, 0.3, 0.4],
        "env_fn": env_scenario_B,
        "param_schedule_fn": None,
    },
    "C_Therapeutic_Recovery": {
        "y0": [0.5, 0.5, 0.4, 0.4, 0.3, 0.4],
        "env_fn": env_scenario_C,
        "param_schedule_fn": scenario_C_params,
    },
}


# ============================================================
# Run all scenarios
# ============================================================
results = {}
for name, cfg in SCENARIOS.items():
    results[name] = run_simulation(
        scenario_name=name,
        y0=cfg["y0"],
        env_fn=cfg["env_fn"],
        params=DEFAULT_PARAMS,
        param_schedule_fn=cfg["param_schedule_fn"],
        t_span=(0.0, 100.0),
        n_eval=1000,
        method="RK45",
    )


# ============================================================
# Plotting
# ============================================================
state_labels = ["T", "D", "S", "E", "R", "M"]
pretty_names = {
    "T": "Tolerance",
    "D": "Danger signal",
    "S": "Cellular stress",
    "E": "Energy reserve",
    "R": "Recovery capacity",
    "M": "Microbiota diversity",
}

for scenario_name, sol in results.items():
    plt.figure(figsize=(10, 6))
    for i, label in enumerate(state_labels):
        plt.plot(sol.t, sol.y[i], label=f"{label} ({pretty_names[label]})")
    plt.title(scenario_name)
    plt.xlabel("Time")
    plt.ylabel("Normalized state")
    plt.ylim(0.0, 1.05)
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()

plt.show()


# ============================================================
# Optional: print final values
# ============================================================
for scenario_name, sol in results.items():
    print(f"\n=== {scenario_name} final state ===")
    final_vals = sol.y[:, -1]
    for label, val in zip(state_labels, final_vals):
        print(f"{label} = {val:.4f}")
