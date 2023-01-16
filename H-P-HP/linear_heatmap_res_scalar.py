import pandas as pd
import numpy as np
import runner
import multiprocessing as mp
import matplotlib.pyplot as plt
from datetime import datetime

def make_linspace(val_min, val_max, resolution):
    
    vals = [val_min + (val_max - val_min)*(i/resolution) for i in range(resolution+1)]

    return vals

def sweep_parameters(resolution, num_cpu, beta_max, beta_lin, alpha_max, alpha_min, sigma_max, b, q, d, rho, etas, gamma, lams, c1, c2, beta_scalar, hyper, seeds, alpha_init, sigma_init, evo_step_count, host_density, para_density, hyper_density, sol):
    
    processes = []
    count = 0
    for i in range(resolution+1):
        for j in range(resolution+1):
            if count < num_cpu:
                p = mp.Process(target = sol.alpha_ad_dyn_v4, args = (beta_max, beta_lin, alpha_max, alpha_min, sigma_max, b, q, d, rho, etas[i], gamma, lams[j], c1, c2, beta_scalar, hyper, seeds[i][j], alpha_init, sigma_init, evo_step_count, host_density, para_density, hyper_density))
                p.start()
                processes.append(p)
                count += 1
            else:
                for proc in processes:
                    proc.join()
                p = mp.Process(target = sol.alpha_ad_dyn_v4, args = (beta_max, beta_lin, alpha_max, alpha_min, sigma_max, b, q, d, rho, etas[i], gamma, lams[j], c1, c2, beta_scalar, hyper, seeds[i][j], alpha_init, sigma_init, evo_step_count, host_density, para_density, hyper_density))
                p.start()
                processes  = [p]
                count = 1    

    for proc in processes:
        proc.join()
    return

def check_branching_flag(seed):
    branching_flag = False
    df = pd.read_csv(f"../data/alpha_evo/data_set{seed}.csv")
    alpha_coords = df["alpha_val"].values
    evo_steps = df["Evolutionary_step"].unique()
    for step in evo_steps[::10]:
        dft = df[df["Evolutionary_step"]==step]
        alpha_coords = sorted(dft["Trait_index_1"].values)
        prev = alpha_coords[0]
        for val in alpha_coords[1:]:
            if (val - prev) > 1:
                branching_flag = True
                break
            prev = val
        if branching_flag:
            break
            
    dft = df[df["Evolutionary_step"]==evo_steps[-1]]
    hyperparasite_flag = 1 - dft["hyperparasites_present"].iloc[0]
    
    return branching_flag, hyperparasite_flag

def create_heatmap_matrice(resolution, res_scalar, num_cpu, max_depth, beta_max, beta_lin, alpha_max, alpha_min, sigma_max, b, q, d, rho, etas, gamma, lams, c1, c2, beta_scalar, hyper, seeds, alpha_init, sigma_init, evo_step_count, host_density, para_density, hyper_density, sol):
    
    print("Doing sweep")
    sweep_parameters(resolution, num_cpu, beta_max, beta_lin, alpha_max, alpha_min, sigma_max, b, q, d, rho, etas, gamma, lams, c1, c2, beta_scalar, hyper, seeds, alpha_init, sigma_init, evo_step_count, host_density, para_density, hyper_density, sol)
    
    heatmap_matrice_temp = []
    print("Making matrice")
    for i, row in enumerate(seeds):
        temp_row = []
        for j, seed in enumerate(row):
            branching_flag, hyperparasite_flag = check_branching_flag(seed)
            temp_row.append(branching_flag*1)
#             temp_row.append(branching_flag*2 + hyperparasite_flag*1)
        heatmap_matrice_temp.append(temp_row)
        
    print(len(heatmap_matrice_temp), len(heatmap_matrice_temp[0]))
#     for row in heatmap_matrice_temp:
#         print(row)
    plt.imshow(heatmap_matrice_temp, origin = "lower")
    plt.savefig("../supplementary_figures/temp_fig.png")
    plt.close()
    
    check_inds = []
    for i in range(resolution):
        for j in range(resolution):
            check_inds.append([i,j])
    
    for depth in range(1,max_depth+1):
        print(depth)
        iis = [0, 1]
        jjs = [0, 1]
        new_etas = []
        new_lams = []
        refine_inds = []
        
        
        for ind in check_inds:
            val = heatmap_matrice_temp[ind[0]][ind[1]]
            flag = False
            for ii in iis:
                for jj in jjs:
                    if heatmap_matrice_temp[ind[0] + ii][ind[1] + jj] != val:
                        flag = True
            if flag:
                    new_etas.append([etas[ind[0]], etas[ind[0]+1]])
                    new_lams.append([lams[ind[1]], lams[ind[1]+1]])
                    refine_inds.append([ind[0], ind[1]])       
                    
#         for ind in refine_inds:
#             for ii in iis:
#                 for jj in jjs:
#                     print(heatmap_matrice_temp[ind[0]  + ii][ind[1] + jj])

        heatmap_matrice_refined = []
        for i, row in enumerate(heatmap_matrice_temp[:-1]):
            for ii in range(res_scalar+1):
                temp_row = []
                for j, val in enumerate(row[:-1]):
                    for jj in range(res_scalar+1):
                        temp_row.append(val)
                temp_row.append(row[-1])
                heatmap_matrice_refined.append(temp_row)
                
#         print(np.shape(heatmap_matrice_temp))
        temp_row = []
        for j, val in enumerate(heatmap_matrice_temp[-1][:-1]):
            for jj in range(res_scalar+1):
                temp_row.append(val)
#         print(len(temp_row))
        temp_row.append(heatmap_matrice_refined[-1][-1])
#         print(len(temp_row))
        heatmap_matrice_refined.append(temp_row)
        
        print(len(heatmap_matrice_refined), len(heatmap_matrice_refined[0]))
        
        for k, ind in enumerate(refine_inds):
            print(f"Done {k+1} of {len(refine_inds)} runs       ", end = "\r")
            eta_pair = new_etas[k]
            lam_pair = new_lams[k]
            temp_etas = make_linspace(eta_pair[0], eta_pair[1], res_scalar)
            temp_lams = make_linspace(lam_pair[0], lam_pair[1], res_scalar)

            sweep_parameters(res_scalar, num_cpu, beta_max, beta_lin, alpha_max, alpha_min, sigma_max, b, q, d, rho, temp_etas, gamma, temp_lams, c1, c2, beta_scalar, hyper, seeds, alpha_init, sigma_init, evo_step_count, host_density, para_density, hyper_density, sol)
            for i in range(res_scalar + 1):
                for j in range(res_scalar + 1):
                    seed = seeds[i][j]
                    branching_flag, hyperparasite_flag = check_branching_flag(seed)
                    heatmap_matrice_refined[ind[0]*(res_scalar + 1) + i][ind[1]*(res_scalar + 1) + j] = 1*branching_flag
                    
#             for i, row in enumerate(seeds):
#                 for j, seed in enumerate(row):
#                     branching_flag, hyperparasite_flag = check_branching_flag(seed)
#                     heatmap_matrice_refined[ind[0]*res_scalar + (i + 1)][ind[1]*res_scalar + (j + 1)] = 1*branching_flag
#                     heatmap_matrice_refined[ind[0]*resolution + (i + 1)][ind[1]*resolution + (j + 1)] = 1*branching_flag + 2*hyperparasite_flag
        heatmap_matrice_temp = heatmap_matrice_refined
        etas = make_linspace(etas[0], etas[-1], len(heatmap_matrice_temp))
        lams = make_linspace(lams[0], lams[-1], len(heatmap_matrice_temp))
        
        check_inds = []
        for ind in refine_inds:
            for i in range(res_scalar + 1):
                for j in range(res_scalar + 1):
                    check_inds.append([ind[0]*res_scalar + i, ind[1]*res_scalar + j])
        
    heatmap_matrice = heatmap_matrice_temp
#     print(heatmap_matrice)
    return heatmap_matrice


t1 = datetime.now()
num_cpu = mp.cpu_count()

sol = runner.PySolver()


b = 2.0
q = 0.1
d = 0.1
gamma = 0.1
sigma = 0.4
sigma_max = sigma
c1 = 1.0

alpha_init = 10
sigma_init = 100
hyper = 1.0
beta_scalar = 1.0
rhos = [0.25, 0.5, 0.75]
lam_min = 0.0
lam_max = 5.0

eta_min = 0.0
eta_max = 2.0
seed_base = 100
evo_step_count = 2000
beta_max = 5.0
alpha_min = 0.0
alpha_max = 5.0
seeds = []

resolution = 30
res_scalar = 10

for k in range(len(rhos)):
    tot_seeds = []
    for i in range(resolution+1):
        temp_seeds = []
        for j in range(resolution+1):
            temp_seeds.append(seed_base + k*((resolution+1)**2) + i*(resolution+1) + j)
        tot_seeds.append(temp_seeds)
    seeds.append(tot_seeds)
    
max_depth = 1
lams = make_linspace(lam_min, lam_max, resolution)
etas = make_linspace(eta_min, eta_max, resolution)

processes = []

# beta_max = 0.5026133332531301
# beta_lin = 0.9157845368366595
# c2 = -68.10004091076856

beta_max = 0.5
beta_lin = 1.0
c2 = -60

seed_for_baseline = 10
rho_baseline = eta_baseline = lam_baseline = 0.0
sol.alpha_ad_dyn_v4(beta_max, beta_lin, alpha_max, alpha_min, sigma_max, b, q, d, rho_baseline, eta_baseline, gamma, lam_baseline, c1, c2, beta_scalar, 0.0, seed_for_baseline, alpha_init, sigma_init)

df_baseline = pd.read_csv(f"../data/alpha_evo/data_set{seed_for_baseline}.csv")
evo_steps = df_baseline["Evolutionary_step"].values

dft = df_baseline[df_baseline["Evolutionary_step"]==evo_steps[-1]]
pops = []
for val in dft["Trait_index_1"].values:
    df_pop = dft[dft["Trait_index_1"]==val][["Density_of_parasite"]]
    pops.append((df_pop["Density_of_parasite"].iloc[0])/(sum(dft["Density_of_parasite"].values)))
    
alpha_weights = []
for i, val in enumerate(pops):
    alpha_weights.append(val*dft["Trait_index_1"].values[i])
alpha_mean = sum(alpha_weights)
# print(alpha_mean)
alpha_init = round(alpha_mean)
# flop

host_density = dft["Density_of_hyperparasite"].iloc[0]
para_density = sum(dft["Density_of_parasite"].values)
hyper_density = 0.1
for k, rho in enumerate(rhos):
    branching_mat = create_heatmap_matrice(resolution, res_scalar, num_cpu, max_depth, beta_max, beta_lin, alpha_max, alpha_min, sigma_max, b, q, d, rho, etas, gamma, lams, c1, c2, beta_scalar, hyper, seeds[k], alpha_init, sigma_init, evo_step_count, host_density, para_density, hyper_density, sol)

    etas_refined = make_linspace(eta_min, eta_max, len(branching_mat) - 1)
    lams_refined = make_linspace(lam_min, lam_max, len(branching_mat) - 1)
#         dists = []
#         for row in branching_mat:
#             dists.append(len(row))
#         print(dists)
#         print(len(branching_mat))
#         print(np.shape(branching_mat))
#         print(np.shape(etas_refined))
#         print(np.shape(lams_refined))
    fig, ax = plt.subplots(1, 1, figsize = (resolution, resolution))
    plt.imshow(branching_mat, origin = "lower")

    ticks_temp = []
    for i in range(len(branching_mat)):
        ticks_temp.append(i)

    ticks = ticks_temp[::resolution]
    ticks.append(ticks_temp[-1])
    lam_ticks = [round(lams_refined[i], 3) for i in ticks]
    eta_ticks = [round(etas_refined[i], 3) for i in ticks]
    plt.yticks(ticks = ticks, labels = eta_ticks, fontsize = max(0.4*resolution, 24))
    plt.xticks(ticks = ticks, labels = lam_ticks, fontsize = max(0.4*resolution, 24))
    plt.title(fr"Probability of cotransmission, $\rho$ = {rho}", fontsize = max(0.5*resolution, 30))
    plt.ylabel(r"$\eta$ value", fontsize = max(0.5*resolution, 30))
    plt.xlabel(r"$\lambda$ value", fontsize = max(0.5*resolution, 30))
    plt.savefig(f"../supplementary_figures/branching_heatmap_{k}_linear_test_comparison.png", bbox_inches = "tight")
    plt.close()
    print(f"The code took {datetime.now() - t1} seconds to run")
    array_for_df = []
    for i, row in enumerate(branching_mat):
        temp_row = [round(etas_refined[i], 3)]
        for val in row:
            temp_row.append(val)
        array_for_df.append(temp_row)
    df = pd.DataFrame(array_for_df)
    df = df.rename(columns = {0: "eta_val"})
    for i, val in enumerate(lams_refined):
        df = df.rename(columns = {i+1: str(round(val,3))})

    df.to_csv(f"../data/branching_mat_{k}_lin_test_comparison.csv", index = False)
            
del sol