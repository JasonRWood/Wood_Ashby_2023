import sys, os, shutil
import runner
import numpy as np
from math import sqrt, exp, dist, floor, ceil, isnan
import plotly.graph_objects as go
from time import sleep
import matplotlib.pyplot as plt
import seaborn as sn
import imageio
import multiprocessing as mp
import pandas as pd

sol = runner.PySolver()

lams = [0.5, 1.0, 2.0]
etas = [0.5, 1.0, 2.0]

b = 2.0
q = 0.1
d = 0.5
gamma = 0.5
hyper = 1.0
c1 = 0.75
c2 = 2.25
rho = 0.5

b_base = 1
q_base = 1

alpha_max = 5
sigma_max = 4
beta_max = 5
c2_min = -2.0
c2_max = 2.0


sigma_min = 0.0
sigma_max_range = 4.0

sigma = 4.0

iter_max = 5
res_scalar = 4
depth_max = iter_max
resolution = 100
param_res = 400

alphas = [alpha_max*i/(param_res) for i in range (param_res + 1)]
betas = [beta_max*i/(param_res) for i in range (param_res + 1)]

seed = 10000

panel_labels = [["(A)", "(B)", "(C)"],["(D)", "(E)", "(F)"], ["(G)","(H)","(I)"]]

fig = plt.figure(figsize = (100, 100))
gs = fig.add_gridspec(3,3)
ax = []

for i in range(3):
    temp_ax = []
    for j in range(3):
        temp_ax.append(fig.add_subplot(gs[i,j]))
    ax.append(temp_ax)
    
for eta_tracker, eta in enumerate(etas):
    for lam_tracker, lam in enumerate(lams):
        current_mat = []
        for beta_val in betas:
            temp_row = []
            for alpha in alphas:
                y = sol.eco_steady_state(
                    beta_val, alpha, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, 0.0, seed
                )
                y = sol.eco_steady_state(
                    beta_val, alpha, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed, y[0], y[1], 1e-3
                )

                S = y[0]
                I = y[1]
                H = y[2]
                
                if S > 1e-3 and I < 1e-3 and H < 1e-3:
                    score = 1
                elif S > 1e-3 and I > 1e-3 and H < 1e-3:
                    score = 2
                elif S > 1e-3 and I > 1e-3 and H > 1e-3:
                    score =3
                else:
                    print(S,I,H)
                    flop
                    
                temp_row.append(score)
            current_mat.append(temp_row)
        
        ax[eta_tracker][lam_tracker].imshow(current_mat, origin = "lower", cmap = "tab10")
        
        ticks_temp = []
        for k in range(len(current_mat)):
            ticks_temp.append(k)

        alpha_ticks = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        beta_ticks = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        ticks = [i for i in ticks_temp if alphas[i] in alpha_ticks]
        
        ax[eta_tracker][lam_tracker].set_yticks(ticks = ticks, labels = beta_ticks, fontsize = 2*resolution/3)
        ax[eta_tracker][lam_tracker].set_xticks(ticks = ticks, labels = alpha_ticks, fontsize = 2*resolution/3, rotation=45)
        ax[eta_tracker][lam_tracker].set_title(fr"$\eta$ = {eta} and $\lambda$ = {lam}", fontsize = 2*resolution/3)
        if lam_tracker == 0 and eta_tracker == 1:
            ax[eta_tracker][lam_tracker].set_ylabel(r"Parasite infectivity, $\beta$", fontsize = 2*resolution/3)
        if lam_tracker == 1 and eta_tracker == 2:
            ax[eta_tracker][lam_tracker].set_xlabel(r"Parasite virulence , $\alpha$", fontsize = 2*resolution/3)
        ax[eta_tracker][lam_tracker].text(0.05,0.9,panel_labels[eta_tracker][lam_tracker], transform=ax[eta_tracker][lam_tracker].transAxes, fontsize = 300/3)    
        
plt.savefig(f"../supplementary_figures/ecological_heatmaps.pdf", bbox_inches = "tight")
plt.savefig(f"../supplementary_figures/ecological_heatmaps.png", bbox_inches = "tight")
plt.close()