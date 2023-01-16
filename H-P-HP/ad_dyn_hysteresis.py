import sys, os, shutil
import runner
import numpy as np

import matplotlib.pyplot as plt


import pandas as pd

rho = 0.9

lam = 2.0

b = 2.0
q = 0.1
d = 0.5
gamma = 0.5

seed = 100

sol = runner.PySolver()

etas = [0.5, 0.2, 0.5, 1.0, 0.5]
    
plotting_alphas = []
plotting_evos = []
alpha_init = 20
sigma_init = 100

beta_max = 5.0
alpha_max = 5.0
sigma_max = 4.0

hyper = 1.0
c1 = 0.75
c2 = 2.25

sol.alpha_ad_dyn(beta_max, alpha_max, sigma_max, b, q, d, rho, 0.0, gamma, 0.0, c1, c2, 0.0, seed, alpha_init, sigma_init, H_density = 0.0)
df = pd.read_csv(f"../data/alpha_evo/data_set{seed}.csv")

evo_steps = df["Evolutionary_step"].values
dft = df[df["Evolutionary_step"]==evo_steps[-1]]
host_density = dft["Density_of_hyperparasite"].iloc[0]
para_density = sum(dft["Density_of_parasite"].values)
hyper_density = 0.1

for eta in etas:
    
    sol.alpha_ad_dyn(beta_max, alpha_max, sigma_max, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed, alpha_init, sigma_init, S_density = host_density, I_density = para_density, H_density = hyper_density)
    df = pd.read_csv(f"../data/alpha_evo/data_set{seed}.csv")

    seed += 1
    
    alpha_coords = df["alpha_val"].values
    evo_steps = df["Evolutionary_step"].values
    
    plotting_alphas.append(alpha_coords)
    plotting_evos.append(evo_steps)
    
    dft = df[df["Evolutionary_step"]==evo_steps[-1]]
    para_densities = dft["Density_of_parasite"].values
    hyper_densities = dft["Density_of_hyperparasite"].values
    trait_inds = dft["Trait_index_1"].values
    alpha_init = 0
    max_density = 0
    for i, den in enumerate(para_densities):
        if (den + hyper_densities[i]) > max_density:
            alpha_init = trait_inds[i]
    
    host_density = dft["Density_of_hyperparasite"].iloc[0]
    para_density = sum(dft["Density_of_parasite"].values)
    hyper_density = sum(dft["Density_of_hyperparasite"].values)
    
del sol


up_shift = []
for i in range(len(etas)):
    up_shift.append(i*1001)
    
colours = {0.5:"black", 0.2:"indigo", 1.0:"red"}
fig, ax = plt.subplots()
    
tagged_labels = []
for i,eta in enumerate(etas):
    
    alphas = plotting_alphas[i]
    evos_temp = plotting_evos[i]
    
    evo_steps = []
    for step in evos_temp:
        evo_steps.append(step + up_shift[i])
        
    if eta not in tagged_labels:
        plt.scatter(alphas, evo_steps, marker = ".", label = fr"$\eta$ = {eta}", c = colours[eta])
        tagged_labels.append(eta)
    else:
        plt.scatter(alphas, evo_steps, marker = ".", c = colours[eta])
        
    
plt.xlabel(r"Parasite virulence, $\alpha$")
plt.xlim([0,3])
plt.ylabel("Evolutionary Time")
ax.legend(loc='center left', bbox_to_anchor=(0.6, 0.85), fontsize = 14)
plt.savefig("../supplementary_figures/ad_dyn_figure.png", bbox_inches='tight')
plt.savefig("../supplementary_figures/ad_dyn_figure.pdf", bbox_inches='tight')
plt.close()