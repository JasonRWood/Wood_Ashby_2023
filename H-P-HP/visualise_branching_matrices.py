import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

rhos = [0.25, 0.5, 0.75]

panel_labels = ["(A)", "(B)", "(C)"]
fig = plt.figure(figsize = (100, 100))
gs = fig.add_gridspec(1,3)
ax = []

for i in range(1):
    temp_ax = []
    for j in range(3):
        temp_ax.append(fig.add_subplot(gs[i,j]))
    ax.append(temp_ax)
    
for i in range(3):
    rho = rhos[i]
    
    df = pd.read_csv(f"../data/branching_mat_{i}_lin_3.csv")
    etas = df["eta_val"].values
    lams = [float(val) for val in df.columns[1:]]
    value_matrix = df[df.columns[1:]].values
    print(value_matrix)
    print(len(etas))
    print(len(lams))
    
    resolution = 100
    
    ax[0][i].imshow(value_matrix, origin = "lower", cmap = "Set3")

    ticks_temp = []
    for k in range(len(value_matrix)):
        ticks_temp.append(k)

#     ticks = ticks_temp[::int(2*resolution)]
#     ticks.append(ticks_temp[-1])
    
    eta_ticks = [0.0, 0.4, 0.8, 1.2, 1.6, 2.0]
    lam_ticks = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    ticks = [i for i in ticks_temp if lams[i] in lam_ticks]
#     lam_ticks = [round(lams[i], 3) for i in ticks]
#     eta_ticks = [round(etas[i], 3) for i in ticks]
    ax[0][i].set_yticks(ticks = ticks, labels = eta_ticks, fontsize = 2*resolution/3)
    ax[0][i].set_xticks(ticks = ticks, labels = lam_ticks, fontsize = 2*resolution/3, rotation=45)
    ax[0][i].set_title(fr"Probability of cotransmission, $\rho$ = {rho}", fontsize = 2*resolution/3)
    if i == 0:
        ax[0][i].set_ylabel(r"Hyperparasite infectivity modifier, $\eta$", fontsize = 2*resolution/3)
    if i == 1:
        ax[0][i].set_xlabel(r"Hyperparasite virulence modifier, $\lambda$", fontsize = 2*resolution/3)
    ax[0][i].text(0.05,0.9,panel_labels[i], transform=ax[0][i].transAxes, fontsize = 300/3)    
plt.savefig(f"../supplementary_figures/branching_heatmaps.pdf", bbox_inches = "tight") 
plt.savefig(f"../supplementary_figures/branching_heatmaps.png", bbox_inches = "tight")
plt.close()