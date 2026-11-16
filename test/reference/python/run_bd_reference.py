# Reference runs of the full host-parasite IPM of Wilber et al. (2022) for test/test_chytrid.jl and
# vignette 24.
#
# This script runs the authors' own code (full_model.py and its import, reduced_model.py, from
# https://github.com/mqwilber/bd_seasonality, commit 4392677, Zenodo doi:10.5281/zenodo.6407088), which is
# not shipped with this package (it carries no licence). Clone that repository, unzip its results/
# pickled_results.zip, and run, with numpy 1.23, scipy, pandas and pyyaml (the code uses np.int and np.bool,
# removed in numpy 1.24):
#
#   python run_bd_reference.py <bd_seasonality>/code <output directory>
#
# The parameters are built as the authors' ccs_simulation.py builds them for its deterministic runs: the
# medians of the MCMC chains of the fitted reduced IPM (K = 4, omega = 1, intensity-dependent loss of
# infection, the model of each location used there), mapped to the full IPM with constant infected
# survival and load-dependent loss, on 300 cells from -20 to 25 ln zoospores, a step of 7 days, a year of
# 364 days and a sinusoidal temperature between the location's minimum and maximum. That mapping is
# re-implemented here (it is a dictionary of renamed values), as is the loading of the medians.
#
# reduced_model.py imports reduced_mcmc.py (arviz, seaborn) and matplotlib, none of which the full model
# uses; empty stand-ins are put in their place.
import sys
import types
import pickle
import csv
import os

code_dir, out_dir = sys.argv[1], sys.argv[2]
for name in ("reduced_mcmc", "matplotlib", "matplotlib.pyplot"):
    sys.modules[name] = types.ModuleType(name)
sys.path.insert(0, code_dir)

import numpy as np          # noqa: E402
import yaml                 # noqa: E402
import full_model as full   # noqa: E402


class _NumpyOnly(pickle.Unpickler):
    """The fitted results hold numpy arrays in dictionaries; refuse anything else."""
    allowed = {("numpy.core.multiarray", "_reconstruct"), ("numpy", "ndarray"), ("numpy", "dtype"),
               ("numpy.core.multiarray", "scalar")}

    def find_class(self, module, name):
        if (module, name) in self.allowed:
            return super().find_class(module, name)
        raise pickle.UnpicklingError("refusing {0}.{1}".format(module, name))


# (location, species, (min, max temperature), model) as in ccs_simulation.py; VT as in model_analysis.ipynb
SITES = [("PA", "LIPI", (4, 25), "a_temp"),
         ("LA", "LISP", (4, 30), "both_model"),
         ("TN", "LISP", (4, 27), "a_temp"),
         ("VT", "LIPI", (4, 25), "a_temp")]
K, OMEGA = 4, 1
IPM = {"min_size": -20, "max_size": 25, "bins": 300, "time_step": 7}
YEAR = 364
AREA = 500.0          # m^2; ccs_simulation.py uses 100 to 1000
DENSITY = 0.14        # post-metamorphic frogs per m^2 at the start, as there
STEPS = 530           # as there: about ten years
SAVED = [0, 1, 26, 52, 104, 260, 529]


def reduced_parameters(dod, spp, model):
    path = os.path.join(code_dir, "..", "results", "pickled_results",
                        "{0}_{3}_{1}_allsites_ext_beta_prior_K={2}_loss_load_omega={4}.pkl".format(dod, spp, K, model, OMEGA))
    with open(path, "rb") as f:
        res = _NumpyOnly(f).load()
    chains, adapt = res["chains"], res["adapt_param"]
    base = dict(res["base_params"][spp])
    draws = np.array([res["mcmc_results"][i]["params"][spp][:, 5 * adapt:].T for i in range(chains)])
    medians = np.percentile(draws, 50, axis=(0, 1))
    for i, p in enumerate(res["parameters"]):
        base[p] = medians[i]
    t = yaml.safe_load(open(os.path.join(code_dir, "model_params", "{0}_params".format(dod),
                                         "spp_params_{0}.yml".format(spp))))
    for hp in ["hibernation", "hibernation_aquatic", "hibernation_start", "hibernation_end", "sigma_full_temp"]:
        base[hp] = t["params"][hp]
    base["base_temp"] = 15
    base["sigmaG"] = base["sigma_full"]
    base["sigma0"] = base["sigma_full"]
    base["omega"] = -1000000
    return base, list(res["parameters"])


def full_parameters(r):
    """The mapping of ccs_simulation.get_full_ipm_params, constant_surv=True, constant_loss=False."""
    return {"s0": r["s0"], "growth_fxn_inter": r["a"], "growth_fxn_temp": r["a_temp"], "growth_fxn_slope": r["b"],
            "growth_fxn_sigma": r["sigmaG"], "loss_fxn_inter": r["mu_l"], "loss_fxn_slope": r["sigma_l"],
            "init_inf_fxn_inter": r["mu0"], "init_inf_fxn_sigma": r["sigma0"], "surv_fxn_inter": r["mu_s"],
            "surv_fxn_slope": r["sigma_s"], "trans_fxn_zpool": r["trans_beta"], "trans_fxn_temp": r["trans_beta_temp"],
            "max_load": 1000, "min_load": -1000, "shedding_prop": r["lam"], "fec": r["r"], "K": r["K"],
            "constant_surv": True, "constant_loss": False, "sI": r["sI"], "lI": 1 / (1 + np.exp(-r["lI"])),
            "larval_period": r["larval_period"], "aquatic": r["aquatic"], "breeding_start": r["breeding_start"],
            "breeding_end": r["breeding_end"], "hibernation": r["hibernation"],
            "hibernation_aquatic": r["hibernation_aquatic"], "hibernation_start": r["hibernation_start"],
            "hibernation_end": r["hibernation_end"], "base_temp": r["base_temp"]}


os.makedirs(out_dir, exist_ok=True)
zsurv = full.zsurv_fxn_temperature(IPM["time_step"])
# the fitted constants of the zoospore survival curve, read back from the function's closure
values = {n: c.cell_contents for n, c in zip(zsurv.__code__.co_freevars, zsurv.__closure__)}

with open(os.path.join(out_dir, "chytrid_zoospore_survival.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["ld50", "s", "maxsurv"])
    w.writerow([repr(float(values["ld50"])), repr(float(values["s"])), repr(float(values["maxsurv"]))])

prow, trow, srow = [], [], []
keys = ["s0", "sI", "growth_fxn_inter", "growth_fxn_temp", "growth_fxn_slope", "growth_fxn_sigma", "loss_fxn_inter",
        "loss_fxn_slope", "trans_fxn_zpool", "trans_fxn_temp", "shedding_prop", "fec", "K", "larval_period",
        "aquatic", "breeding_start", "breeding_end", "hibernation", "hibernation_aquatic", "hibernation_start",
        "hibernation_end", "base_temp"]
for dod, spp, (tmin, tmax), model in SITES:
    red, fitted = reduced_parameters(dod, spp, model)
    pf = full_parameters(red)
    prow.append([dod, spp, model, " ".join(fitted)] + [repr(float(pf[k])) for k in keys])
    init = np.r_[0, DENSITY * AREA, np.zeros(IPM["bins"]), 10000.0]
    comm = full.Community(dod, {"time": 0, "species": {spp: pf}, "density": {spp: init}}, IPM, zsurv,
                          area=AREA, base_temp=pf["base_temp"], init_temp=tmin, min_temp=tmin, max_temp=tmax,
                          year_length=YEAR)
    sp = comm.species[0]
    sim = np.empty((IPM["bins"] + 3, STEPS))
    sim[:, 0] = sp.density_vect
    temperature = np.empty(STEPS)
    temperature[0] = comm.temperature
    for step in range(1, STEPS):
        comm.update_deterministic()
        sim[:, step] = sp.density_vect
        temperature[step] = comm.temperature
    y = sp.y
    for k in range(STEPS):
        v = sim[:, k]
        infected = v[2:-1]
        ni = infected.sum()
        mean = float((infected * y).sum() / ni) if ni > 0 else float("nan")
        trow.append([dod, k, repr(float(temperature[k])), repr(float(v[0])), repr(float(v[1])), repr(float(ni)),
                     repr(float(v[-1])), repr(mean)])
    for k in SAVED:
        srow.append([dod, k] + [repr(float(x)) for x in sim[:, k]])

with open(os.path.join(out_dir, "chytrid_parameters.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["location", "species", "model", "fitted"] + keys)
    w.writerows(prow)
with open(os.path.join(out_dir, "chytrid_runs.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["location", "step", "temperature", "tadpoles", "susceptible", "infected", "zoospores", "mean_load"])
    w.writerows(trow)
with open(os.path.join(out_dir, "chytrid_states.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["location", "step", "values: tadpoles, susceptible, 300 infected cells, zoospores"])
    w.writerows(srow)
print("done")
