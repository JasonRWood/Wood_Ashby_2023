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

sys.path.append("../src/Models/")

import singular_strategy as ss

sys.path.append("../src/zero_finder/")

import zero_finder_utilities as zf

sys.path.append("../src/figures/")

import figure_utilities as ff

from datetime import datetime

with open("../data/words.txt", encoding='utf8') as f:
    lines = f.readlines()
    
lines_cleaned = []

for line in lines:
    if line[0] != "#" and "www" not in line and line[0] != " ":
        lines_cleaned.append(line[:-1])

def calculate_alpha_gradient_vector(sol, resolution, beta_max, alpha_max, alpha_min_res, alpha_max_res, sigma_max, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed):
    
    alpha_grad_vec = []
    alpha_val = []
    
    for i in range(resolution+1):

        alpha = alpha_min_res + (alpha_max_res - alpha_min_res) * (i / resolution)
        alpha_val.append(alpha)

        beta = zf.beta_func(beta_max, alpha, alpha_max, sigma, sigma_max, c1, c2)
        
        dbetadalpha = zf.dbetadalpha_func(
            beta_max, alpha, alpha_max, sigma, sigma_max, c1, c2
        )

        y = sol.eco_steady_state(
            beta, alpha, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed
        )

        S = y[0]
        I = y[1]
        H = y[2]

        if (I == 0 and H == 0):
            alpha_grad = 1*(alpha <= alpha_max/2) - 1*(alpha> alpha_max/2)
        elif (H == 0):
            alpha_grad = dbetadalpha*S/(d + alpha + gamma) - (beta*S)/((d + alpha + gamma)**2)
        elif (I == 0):
            alpha_grad = dbetadalpha*S/(d + lam*alpha + gamma) - (beta*S)/((d + lam*alpha + gamma)**2)
        else:
            alpha_grad = ss.fitness_gradient_alpha(
                beta, alpha, sigma, dbetadalpha, H, S, eta, lam, gamma, d, rho
            )

        alpha_grad = (alpha_grad > 0)*1 + (alpha_grad < 0)*-1
        alpha_grad_vec.append(alpha_grad)
        
    return alpha_grad_vec, alpha_val

def return_zeros(grad_vec, param_space):
    
    poss_attractors = []
    poss_repellers = []
    for i,val in enumerate(grad_vec[:-1]):
        if val == 1 and grad_vec[i+1] == -1:
            poss_attractors.append([param_space[i], param_space[i+1]])
        elif val == -1 and grad_vec[i+1] == 1:
            poss_repellers.append([param_space[i], param_space[i+1]])
        elif val == 0 and grad_vec[i+1] == -1:
            poss_attractors.append([param_space[i], param_space[i+1]])
            
    return poss_attractors, poss_repellers 
    
def get_sing_strats(grad_vec, param_space, iter_max, res_scalar, sol, resolution, beta_max, alpha_max, alpha_min_res, alpha_max_res, sigma_max, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed):
    
    sing_vals = []
    
    if grad_vec[0] == -1:
        sing_vals.append(param_space[0])
    if grad_vec[-1] == 1:
        sing_val.append(param_space[-1])
       
    param_space_temp = param_space
    poss_zeros = return_zeros(grad_vec,param_space_temp)
    for pair in poss_zeros:
        grad_vec_temp, param_space_temp = calculate_alpha_gradient_vector(sol, res_scalar, beta_max, alpha_max, pair[0], pair[1], sigma_max, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed)
        
    
    
    return sing_vals

def recusive_sing_strat_function(sol, depth, iter_max, resolution, beta_max, alpha_max, alpha_min_res, alpha_max_res, sigma_max, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed):
    
    sing_strats = []
    
    attractors = []
    repellers = []
    
    if depth == iter_max:
        
        alpha_grad_vec, alpha_val = calculate_alpha_gradient_vector(sol, resolution, beta_max, alpha_max, alpha_min_res, alpha_max_res, sigma_max, sigma, b, q, d, rho, eta, gamma, lam,c1, c2, hyper, seed)
        
        poss_attractors, poss_repellers = return_zeros(alpha_grad_vec,alpha_val)
        
        for zeros in poss_attractors:
#             attractors.append((zeros[1] + zeros[0])/2)
            attractors.append(zeros[0])
        for zeros in poss_repellers:
#             repellers.append((zeros[1] + zeros[0])/2)
            repellers.append(zeros[0])
            
    else:
        alpha_grad_vec, alpha_val = calculate_alpha_gradient_vector(sol, resolution, beta_max, alpha_max, alpha_min_res, alpha_max_res, sigma_max, sigma, b, q, d, rho, eta, gamma, lam,c1, c2, hyper, seed)
        
        poss_attractors, poss_repellers = return_zeros(alpha_grad_vec,alpha_val)
        
        if depth == 0:
            if alpha_grad_vec[0] == -1:
                attractors.append(alpha_val[0])
            if alpha_grad_vec[-1] == 1:
                attractors.append(alpha_val[-1])
                
        for pair in poss_attractors:
            temp_attractors, temp_repellers = recusive_sing_strat_function(sol, depth+1, iter_max, resolution, beta_max, alpha_max, pair[0], pair[1], sigma_max, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed)
            
            for val in temp_attractors:
                attractors.append(val)
                
        for pair in poss_repellers:
            temp_attractors, temp_repellers = recusive_sing_strat_function(sol, depth+1, iter_max, resolution, beta_max, alpha_max, pair[0], pair[1], sigma_max, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed)
                
            for val in temp_repellers:
                repellers.append(val)
                
    if depth == 0:
        attractors = zf.clean_up_sing_strats(attractors)
        repellers = zf.clean_up_sing_strats(repellers)
                
    return attractors, repellers

def produce_suplimentary_data(
    sol,
    plotting_sing_strats_hyper,
    plotting_hyper,
    plotting_sing_strats_no_hyper,
    plotting_no_hyper,
    beta_max,
    alpha_max,
    sigma_max,
    sigma,
    c1,
    c2,
    b,
    q,
    d,
    rho,
    gamma,
    lam,
    hyper,
    seed,
):

    death_measure_1 = []

    death_measure_2 = []

    relative_pop_size = []

    relative_infected_proportion = []

    hyperparasite_densities = []
    infected_densities = []
    host_densities = []
    
    for i, alpha in enumerate(plotting_sing_strats_hyper):

        eta = plotting_hyper[i]
        no_hyper_ind = plotting_no_hyper.index(eta)

        alpha_no_hyper = plotting_sing_strats_no_hyper[no_hyper_ind]

        beta_val = zf.beta_func(beta_max, alpha, alpha_max, sigma, sigma_max, c1, c2)
        beta_no_hyper = zf.beta_func(beta_max, alpha_no_hyper, alpha_max, sigma, sigma_max, c1, c2)

        y = sol.eco_steady_state(
            beta_val, alpha, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed
        )

        S = y[0]
        I = y[1]
        H = y[2]
        
        host_densities.append(S)
        infected_densities.append(I)
        hyperparasite_densities.append(H)
        
#         while H == 0 and I == 0 and alpha != alpha_no_hyper:
#             alpha = alpha - 1e-2
#             beta_val = zf.beta_func(beta_max, alpha, alpha_max, sigma, sigma_max, c1, c2)
#             beta_no_hyper = zf.beta_func(beta_max, alpha_no_hyper, alpha_max, sigma, sigma_max, c1, c2)

#             y = sol.eco_steady_state(
#                 beta_val, alpha, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed
#             )

#             S = y[0]
#             I = y[1]
#             H = y[2]
        N = S + I + H 

        y_no_hyper = sol.eco_steady_state(
            beta_no_hyper, alpha_no_hyper, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, 0.0, seed
        )

        S_no_hyper = y_no_hyper[0]
        I_no_hyper = y_no_hyper[1]
        H_no_hyper = y_no_hyper[2]

        N_no_hyper = S_no_hyper + I_no_hyper + H_no_hyper

        deaths_1 = ((d*(S + I + H) + alpha*I + lam*alpha*H)/(S+I+H))/((d*(S_no_hyper + I_no_hyper) + alpha_no_hyper*I_no_hyper)/(S_no_hyper + I_no_hyper))

        if I + H == 0:
            deaths_2 = 0
        else:
            deaths_2 = ((alpha*I + lam*alpha*H)/(I+H))/alpha_no_hyper

        death_measure_1.append(deaths_1)
        death_measure_2.append(deaths_2)

        relative_pop_size.append(N/N_no_hyper)

        if I_no_hyper + H_no_hyper == 0 and I + H != 0:
            relative_infected_proportion.append(10)
        elif I_no_hyper + H_no_hyper == 0 and I + H == 0:
            relative_infected_proportion.append(0)
        else:
            relative_infected_proportion.append(((I + H)/N)/((I_no_hyper + H_no_hyper)/N_no_hyper))
            
    return death_measure_1, death_measure_2, relative_pop_size, relative_infected_proportion, host_densities, infected_densities, hyperparasite_densities

# rhos = [0.0, 0.5, 1.0]
rhos = [0.1, 0.5, 0.9]
# rhos = [0.85]

# bs = [1.5, 2.0, 2.5, 3.0, 4.0]
# qs = [0.1, 0.25, 0.4]
# ds = [0.1, 0.2, 0.5]
# gammas = [0.5, 0.75, 1.0]
# lams = [3, 4, 5, 6]
# lams = [0.05, 0.5, 1, 1.5, 2.5, 5]
lams = [0.5, 1, 2.0]
# lams = [1.0]
b = 2.0
q = 0.1
d = 0.5
gamma = 0.5
hyper = 1.0
c1 = 0.75
c2 = 2.25
b_base = 1
q_base = 1

alpha_max = 5/b_base
sigma_max = 4/q_base
beta_max = 5/q_base
c2_min = -2.0
c2_max = 2.0
# lams = [2.45]

sigma_min = 0.0
sigma_max_range = 4.0
# sigmas = [sigma_min + (sigma_max_range - sigma_min)*(i/10) for i in range(11)]
sigma = 4.0

iter_max = 5
res_scalar = 4
depth_max = iter_max
resolution = 100
param_res = 200
# bs = [1.0]
# qs = [1.0]
# ds = [0.5]

seed = 10000

sol = runner.PySolver()

plotting_dict_attractors = {}
plotting_dict_attractors_no_hyper = {}

plotting_dict_repellers = {}
plotting_dict_repellers_no_hyper = {}

plotting_dict_hyper = {}
plotting_dict_no_hyper = {}
plotting_dict_deaths_1 = {}
plotting_dict_deaths_2 = {}
plotting_dict_pops = {}
plotting_dict_inf = {}
plotting_dict_hosts = {}
plotting_dict_paras = {}
plotting_dict_hypers = {}

labels_base = []
labels = []
lines1 = []
lines2 = []

inds = np.random.randint(len(lines_cleaned),size = 4)

filekey = ""
for ind in inds:
    filekey += lines_cleaned[ind] + "_"

filekey += "alpha_figs"

filekey = filekey.replace("'", "").replace(".", "")

output_folder = "../outputs/single_evo/" + filekey

os.makedirs(output_folder)

with open(f"{output_folder}/parameters.txt", "w") as f:
    f.write("The parameters for the system where:\n")
    f.write(f"b = {b}\n")
    f.write(f"q = {q}\n")
    f.write(f"d = {d}\n")
    for i, rho in enumerate(rhos):
        f.write(f"rho_{i} = {rho}\n")
    for i, lam in enumerate(lams):
        f.write(f"lam_{i} = {lam}\n")
    f.write(f"gamma = {gamma}\n")
    f.write(f"eta = varying\n")
    f.write(f"c1 = {c1}\n")
    f.write(f"c2 = {c2}\n")
    f.write(f"beta_max = {beta_max}\n")
    f.write(f"alpha_max = {alpha_max}\n")
    f.write(f"sigma_max = {sigma_max}\n")
    f.write(f"resolution = {resolution}\n")
    f.write(f"param_res = {param_res}\n")
    f.write(f"iter_max = {iter_max}\n")
    f.write(f"res_scalar = {res_scalar}\n")
    f.write(f"sigma = {sigma}\n")

f.close()

fig1 = plt.figure(figsize = (40, 10))
gs1 = fig1.add_gridspec(1,3)
ax1 = []

fig2 = plt.figure(figsize = (40, 20))
gs2 = fig2.add_gridspec(2,3)
ax2 = []

texts = [["(A)", "(B)", "(C)"],["(D)", "(E)", "(F)"], ["(G)","(H)","(I)"]]

titles = ["Hypovirulence ", "No effect on virulence ", "Hypervirulence "]

for i in range(1):
    temp_ax = []
    for j in range(3):
        temp_ax.append(fig1.add_subplot(gs1[i,j]))
    ax1.append(temp_ax)
    
for i in range(2):
    temp_ax = []
    for j in range(3):
        temp_ax.append(fig2.add_subplot(gs2[i,j]))
    ax2.append(temp_ax)
    
for lam_tracker, lam in enumerate(lams):
    for rho in rhos:


        b_base = 1
        q_base = 1

        alpha_max = 5/b_base
        sigma_max = sigma/q_base
        beta_max = 5/q_base

        sigma_res = [sigma_max*(i/100) for i in range(101)]

        for i, val in enumerate(sigma_res):
            if val >= sigma:
                sigma_init = i
                break

        iter_max = 5
        res_scalar = 4
        depth_max = iter_max
        resolution = 100

        beta_max = beta_max
        alpha_max_param = alpha_max 
        sigma_max_param = sigma_max

        param_res = 200

#                 b = 1.0
#                 q = 1.0

#         d = d_temp/b_base

        eta = 1.0
#                         lam = 1.0

        eta_min = 0.0
        eta_max = 2.0
        etas = [eta_min + (eta_max - eta_min)*i/param_res for i in range(param_res + 1)]

#                             etas = [0.44]
#                             sigma = sigma_max

        
        c2s = [c2_min + (c2_max - c2_min)*(i/param_res) for i in range(param_res+1)]

#                     gamma = 0.5/b_base

        gamma_min = 0.0/b_base
        gamma_max = 1.0/b_base
        gammas = [gamma_min + (gamma_min - gamma_max)*(i/param_res) for i in range(param_res+1)]

        plotting_attractors_hyper = []
        plotting_attractors_no_hyper = []

        plotting_repellers_hyper = []
        plotting_repellers_no_hyper = []

        plotting_etas_attractors_hyper = []
        plotting_etas_attractors_no_hyper = []

        plotting_etas_repellers_hyper = []
        plotting_etas_repellers_no_hyper = []

#                             rho_str = str(rho)
#                             rho_str = rho_str.replace(".", "")
#                             temp_folder_path = f"../outputs/single_evo/temp_folder_{rho_str}/"
#                             os.makedirs(temp_folder_path)

        for eta_t in etas:
            eta = round(eta_t,4)
            print(f"On run {eta} on rho {rho}  ", end = "\r")

            alpha_max_res = alpha_max
            alpha_min_res = 0.0

            alpha_grad_vec, alpha_val = calculate_alpha_gradient_vector(sol, resolution, beta_max, alpha_max, alpha_min_res, alpha_max_res, sigma_max, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed)


            alpha_grad_no_hyper_vec, alpha_val = calculate_alpha_gradient_vector(sol, resolution, beta_max, alpha_max, alpha_min_res, alpha_max_res, sigma_max, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, 0.0, seed)


            attractors_hyper, repellers_hyper = recusive_sing_strat_function(sol, 0, iter_max, resolution, beta_max, alpha_max, alpha_min_res, alpha_max_res, sigma_max, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, hyper, seed)

            attractors_no_hyper, repellers_no_hyper = recusive_sing_strat_function(sol, 0, iter_max, resolution, beta_max, alpha_max, alpha_min_res, alpha_max_res, sigma_max, sigma, b, q, d, rho, eta, gamma, lam, c1, c2, 0.0, seed)

            for val in attractors_hyper:
                plotting_etas_attractors_hyper.append(eta)
                plotting_attractors_hyper.append(val)

            for val in attractors_no_hyper:
                plotting_etas_attractors_no_hyper.append(eta)
                plotting_attractors_no_hyper.append(val)

            for val in repellers_hyper:
                plotting_etas_repellers_hyper.append(eta)
                plotting_repellers_hyper.append(val)

            for val in repellers_no_hyper:
                plotting_etas_repellers_no_hyper.append(eta)
                plotting_repellers_no_hyper.append(val)

            for i, val in enumerate(alpha_val):
                if val >= attractors_no_hyper[0]:
                    break

        plotting_dict_hyper[rho] = {}
        plotting_dict_no_hyper[rho] = {}

        plotting_dict_hyper[rho]["a"] = plotting_etas_attractors_hyper
        plotting_dict_no_hyper[rho]["a"] = plotting_etas_attractors_no_hyper

        plotting_dict_hyper[rho]["r"] = plotting_etas_repellers_hyper
        plotting_dict_no_hyper[rho]["r"] = plotting_etas_repellers_no_hyper

        plotting_dict_attractors[rho] = plotting_attractors_hyper
        plotting_dict_attractors_no_hyper[rho] = plotting_attractors_no_hyper

        plotting_dict_repellers[rho] = plotting_repellers_hyper
        plotting_dict_repellers_no_hyper[rho] = plotting_repellers_no_hyper

        death_measure_1_attractor, death_measure_2_attractor, relative_pop_size_attractor, relative_infected_proportion_attractor, host_density_attractor, para_density_attractor, hyper_density_attractor = produce_suplimentary_data(
            sol,
            plotting_attractors_hyper,
            plotting_etas_attractors_hyper,
            plotting_attractors_no_hyper,
            plotting_etas_attractors_no_hyper,
            beta_max,
            alpha_max,
            sigma_max,
            sigma,
            c1,
            c2,
            b,
            q,
            d,
            rho,
            gamma,
            lam,
            hyper,
            seed,
          )

        death_measure_1_repeller, death_measure_2_repeller, relative_pop_size_repeller, relative_infected_proportion_repeller, host_density_repeller, para_density_repeller, hyper_density_repeller = produce_suplimentary_data(
            sol,
            plotting_repellers_hyper,
            plotting_etas_repellers_hyper,
            plotting_attractors_no_hyper,
            plotting_etas_attractors_no_hyper,
            beta_max,
            alpha_max,
            sigma_max,
            sigma,
            c1,
            c2,
            b,
            q,
            d,
            rho,
            gamma,
            lam,
            hyper,
            seed,
          )

        plotting_dict_deaths_1[rho] = {}
        plotting_dict_deaths_2[rho] = {}
        plotting_dict_pops[rho] = {}
        plotting_dict_inf[rho] = {}

        plotting_dict_deaths_1[rho]["a"] = death_measure_1_attractor
        plotting_dict_deaths_2[rho]["a"] = death_measure_2_attractor
        plotting_dict_pops[rho]["a"] = relative_pop_size_attractor
        plotting_dict_inf[rho]["a"] = relative_infected_proportion_attractor

        plotting_dict_deaths_1[rho]["r"] = death_measure_1_repeller
        plotting_dict_deaths_2[rho]["r"] = death_measure_2_repeller
        plotting_dict_pops[rho]["r"] = relative_pop_size_repeller
        plotting_dict_inf[rho]["r"] = relative_infected_proportion_repeller

        plotting_dict_hosts[rho] = {}
        plotting_dict_paras[rho] = {}
        plotting_dict_hypers[rho] = {}

        plotting_dict_hosts[rho]["a"] = host_density_attractor
        plotting_dict_paras[rho]["a"] = para_density_attractor
        plotting_dict_hypers[rho]["a"] = hyper_density_attractor

        plotting_dict_hosts[rho]["r"] = host_density_repeller
        plotting_dict_paras[rho]["r"] = para_density_repeller
        plotting_dict_hypers[rho]["r"] = hyper_density_repeller


                        

    xlabel = r"Hyperparasite infectivity modifier, $\eta$"
    colours = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:gray"]

    eta_mins = {}
    eta_maxs = {}
    for rho in rhos:
        temp_etas = plotting_dict_hyper[rho]["r"]
        if len(temp_etas) > 0:
            eta_mins[rho] = min(temp_etas)
            eta_maxs[rho] = max(temp_etas)

    alpha_lower_plotting = {}
    alpha_upper_plotting = {}
    alpha_attractors = {}
    eta_attractors = {}
    for rho in rhos:
        temp_alphas = plotting_dict_repellers[rho]
        eta_attractors[rho] = {}
        alpha_attractors[rho] = {}
        eta_attractors[rho]['1'] = []
        eta_attractors[rho]['2'] = []
        alpha_attractors[rho]['1'] = []
        alpha_attractors[rho]['2'] = []
        if len(temp_alphas) > 0:
            alpha_lower_plotting[rho] = min(temp_alphas)
            alpha_upper_plotting[rho] = max(temp_alphas)
            
            flag = True
            paired_list = []
            for i, val in enumerate(plotting_dict_attractors[rho]):
                paired_list.append([val, plotting_dict_hyper[rho]["a"][i]])
            sorted_attractors = sorted(paired_list, key = lambda x: x[0], reverse = True)
            temp_attractors_1 = []
            temp_attractors_2 = []
            temp_etas_1 = []
            temp_etas_2 = []
            for i, val in enumerate(sorted_attractors[:-1]):
                if flag and abs(sorted_attractors[i+1][0] - val[0]) < 0.1:
                    temp_attractors_1.append(val[0])
                    temp_etas_1.append(val[1])
                elif flag:
                    flag = False
                    temp_attractors_1.append(val[0])
                    temp_etas_1.append(val[1])
                else:
                    temp_attractors_2.append(val[0])
                    temp_etas_2.append(val[1])
                    
            temp_attractors_2.append(sorted_attractors[-1][0])
            temp_etas_2.append(sorted_attractors[-1][1])
            
            paired_list_1 = []
            for i, val in enumerate(temp_etas_1):
                paired_list_1.append([val, temp_attractors_1[i]])
             
            paired_list_2 = []
            for i, val in enumerate(temp_etas_2):
                paired_list_2.append([val, temp_attractors_2[i]])
                
            sorted_list_1 = sorted(paired_list_1, key = lambda x: x[0])
            sorted_list_2 = sorted(paired_list_2, key = lambda x: x[0])
            
            temp_attractors_1 = [val[1] for val in sorted_list_1]
            temp_attractors_2 = [val[1] for val in sorted_list_2]
            
            temp_etas_1 = [val[0] for val in sorted_list_1]
            temp_etas_2 = [val[0] for val in sorted_list_2]
            
            alpha_attractors[rho]['1'] = temp_attractors_1
            alpha_attractors[rho]['2'] = temp_attractors_2
            eta_attractors[rho]['1'] = temp_etas_1
            eta_attractors[rho]['2'] = temp_etas_2
            
            plotting_dict_deaths_1[rho]["1"], plotting_dict_deaths_2[rho]["1"], plotting_dict_pops[rho]["1"], plotting_dict_inf[rho]["1"], host_density_attractor, para_density_attractor, hyper_density_attractor = produce_suplimentary_data(
            sol,
            temp_attractors_1,
            temp_etas_1,
            plotting_attractors_no_hyper,
            plotting_etas_attractors_no_hyper,
            beta_max,
            alpha_max,
            sigma_max,
            sigma,
            c1,
            c2,
            b,
            q,
            d,
            rho,
            gamma,
            lam,
            hyper,
            seed,
          )
            
            plotting_dict_deaths_1[rho]["2"], plotting_dict_deaths_2[rho]["2"], plotting_dict_pops[rho]["2"], plotting_dict_inf[rho]["2"], host_density_attractor, para_density_attractor, hyper_density_attractor = produce_suplimentary_data(
            sol,
            temp_attractors_2,
            temp_etas_2,
            plotting_attractors_no_hyper,
            plotting_etas_attractors_no_hyper,
            beta_max,
            alpha_max,
            sigma_max,
            sigma,
            c1,
            c2,
            b,
            q,
            d,
            rho,
            gamma,
            lam,
            hyper,
            seed,
          )
        else:
            alpha_attractors[rho]['1'] = plotting_dict_attractors[rho]
            eta_attractors[rho]['1'] = plotting_dict_hyper[rho]["a"]

            plotting_dict_deaths_1[rho]["1"], plotting_dict_deaths_2[rho]["1"], plotting_dict_pops[rho]["1"], plotting_dict_inf[rho]["1"], host_density_attractor, para_density_attractor, hyper_density_attractor = produce_suplimentary_data(
            sol,
            plotting_dict_attractors[rho],
            plotting_dict_hyper[rho]["a"],
            plotting_attractors_no_hyper,
            plotting_etas_attractors_no_hyper,
            beta_max,
            alpha_max,
            sigma_max,
            sigma,
            c1,
            c2,
            b,
            q,
            d,
            rho,
            gamma,
            lam,
            hyper,
            seed,
          )
            plotting_dict_deaths_1[rho]["2"] = []
            plotting_dict_deaths_2[rho]["2"] = []
            plotting_dict_pops[rho]["2"] = []
            plotting_dict_inf[rho]["2"] = []
            
    for i, rho in enumerate(rhos):
        if lam_tracker == 0:
            labels_base.append(fr"$\rho$ = {rho}")
            labels.append(labels_base[i])
            l1 = ax1[0][lam_tracker].plot(eta_attractors[rho]['1'], alpha_attractors[rho]['1'], c = f"{colours[i]}")
            lines1.append(l1)
            
    for i, rho in enumerate(rhos):
        if lam_tracker == 0:
            l1 = ax2[0][lam_tracker].plot(eta_attractors[rho]['1'], plotting_dict_deaths_2[rho]["1"], c = f"{colours[i]}")
            lines2.append(l1)        
            
    for i, rho in enumerate(rhos):
        if lam_tracker == 0:
#             labels_base.append(fr"$\rho$ = {rho}")
#             labels.append(labels_base[i])
#             l1 = ax[0][lam_tracker].plot(eta_attractors[rho]['1'], alpha_attractors[rho]['1'], c = f"{colours[i]}")
            ax1[0][lam_tracker].plot(eta_attractors[rho]['1'], alpha_attractors[rho]['1'], c = f"{colours[i]}")
#             lines.append(l1)
#             ax[0][lam_tracker].plot(plotting_dict_hyper[rho]["r"], plotting_dict_repellers[rho],  c = f"{colours[i]}", linestyle='dashed')
            ax1[0][lam_tracker].plot(eta_attractors[rho]['2'], alpha_attractors[rho]['2'], c = f"{colours[i]}")
            ax1[0][lam_tracker].scatter(eta_attractors[rho]['1'][-1], alpha_attractors[rho]['1'][-1], edgecolors = f"{colours[i]}", marker = "o", facecolors= f"{colours[i]}", s = 200)
            ax1[0][lam_tracker].scatter(eta_attractors[rho]['1'][0], alpha_attractors[rho]['1'][0], edgecolors = f"{colours[i]}", marker = "o", facecolors="none", s = 200)
            try:
                ax1[0][lam_tracker].scatter(eta_attractors[rho]['2'][-1], alpha_attractors[rho]['2'][-1], edgecolors = f"{colours[i]}", marker = "s", facecolors= f"{colours[i]}", s = 200)
                ax1[0][lam_tracker].scatter(eta_attractors[rho]['2'][0], alpha_attractors[rho]['2'][0], edgecolors = f"{colours[i]}", marker = "s", facecolors="none", s = 200)
            except:
                pass
            if len(plotting_dict_hyper[rho]["r"]) > 0:
                repeller_etas = [eta_attractors[rho]["2"][0]]
                for val in plotting_dict_hyper[rho]["r"]:
                    repeller_etas.append(val)
                repeller_etas.append(eta_attractors[rho]["1"][-1])

                repeller_alphas = [alpha_attractors[rho]['2'][0]]
                for val in plotting_dict_repellers[rho]:
                    repeller_alphas.append(val)

                repeller_alphas.append(alpha_attractors[rho]["1"][-1])
            else:
                repeller_etas = plotting_dict_hyper[rho]["r"]
                repeller_alphas = plotting_dict_repellers[rho]
#             ax[0][lam_tracker].plot(plotting_dict_hyper[rho]["r"], plotting_dict_repellers[rho],  c = f"{colours[i]}", linestyle='dashed')
            ax1[0][lam_tracker].plot(repeller_etas, repeller_alphas,  c = f"{colours[i]}", linestyle='dashed')
        else:
            ax1[0][lam_tracker].plot(eta_attractors[rho]['1'], alpha_attractors[rho]['1'], c = f"{colours[i]}")
            ax1[0][lam_tracker].plot(eta_attractors[rho]['2'], alpha_attractors[rho]['2'], c = f"{colours[i]}")
            ax1[0][lam_tracker].scatter(eta_attractors[rho]['1'][0], alpha_attractors[rho]['1'][0], edgecolors = f"{colours[i]}", marker = "o", facecolors="none", s = 200)
            ax1[0][lam_tracker].scatter(eta_attractors[rho]['1'][-1], alpha_attractors[rho]['1'][-1], edgecolors = f"{colours[i]}", marker = "o", facecolors= f"{colours[i]}", s = 200)
            try:
                ax1[0][lam_tracker].scatter(eta_attractors[rho]['2'][0], alpha_attractors[rho]['2'][0], edgecolors = f"{colours[i]}", marker = "s", facecolors = "none", s = 200)
                ax1[0][lam_tracker].scatter(eta_attractors[rho]['2'][-1], alpha_attractors[rho]['2'][-1], edgecolors = f"{colours[i]}", marker = "s", facecolors= f"{colours[i]}", s = 200)
            except:
                pass
            if len(plotting_dict_hyper[rho]["r"]) > 0:
                repeller_etas = [eta_attractors[rho]["2"][0]]
                for val in plotting_dict_hyper[rho]["r"]:
                    repeller_etas.append(val)
                repeller_etas.append(eta_attractors[rho]["1"][-1])

                repeller_alphas = [alpha_attractors[rho]['2'][0]]
                for val in plotting_dict_repellers[rho]:
                    repeller_alphas.append(val)

                repeller_alphas.append(alpha_attractors[rho]["1"][-1])
            else:
                repeller_etas = plotting_dict_hyper[rho]["r"]
                repeller_alphas = plotting_dict_repellers[rho]
#             ax[0][lam_tracker].plot(plotting_dict_hyper[rho]["r"], plotting_dict_repellers[rho],  c = f"{colours[i]}", linestyle='dashed')
            ax1[0][lam_tracker].plot(repeller_etas, repeller_alphas,  c = f"{colours[i]}", linestyle='dashed')


    ax1[0][lam_tracker].plot(plotting_dict_no_hyper[rho]["a"],plotting_dict_attractors_no_hyper[rho], c = "k")
    ax1[0][lam_tracker].plot(plotting_dict_no_hyper[rho]["r"],plotting_dict_repellers_no_hyper[rho], c = "k")
    ax1[0][lam_tracker].text(0.05,0.95,texts[0][lam_tracker], transform=ax1[0][lam_tracker].transAxes, fontsize = 30)
#     ax[0][lam_tracker].set_xlabel(xlabel, fontsize = 24)
    ax1[0][lam_tracker].set_xlabel("", fontsize = 0)
    ax1[0][lam_tracker].set_ylim([0, alpha_max_res])
    ax1[0][lam_tracker].tick_params(axis='both', which='major', labelsize=34)
    if lam_tracker == 0:
        ax1[0][lam_tracker].set_ylabel(r"Evolved virulence ($\alpha$) value", fontsize = 34)
    ax1[0][lam_tracker].set_title(fr"{titles[lam_tracker]}($\lambda$ = {lam})", fontsize = 34)
    
    if lam_tracker == 1:
        ax1[0][lam_tracker].set_xlabel(xlabel, fontsize = 34)
    else:
        ax1[0][lam_tracker].set_xlabel("", fontsize = 0)
        
    for i, rho in enumerate(rhos):
        ax2[0][lam_tracker].plot(eta_attractors[rho]['1'], plotting_dict_deaths_2[rho]["1"], c = f"{colours[i]}")
        ax2[0][lam_tracker].plot(eta_attractors[rho]['2'], plotting_dict_deaths_2[rho]["2"], c = f"{colours[i]}")
        ax2[0][lam_tracker].scatter(eta_attractors[rho]['1'][0], plotting_dict_deaths_2[rho]['1'][0], edgecolors = f"{colours[i]}", marker = "o", facecolors="none", s = 200)
        ax2[0][lam_tracker].scatter(eta_attractors[rho]['1'][-1], plotting_dict_deaths_2[rho]['1'][-1], edgecolors = f"{colours[i]}", marker = "o", facecolors= f"{colours[i]}", s = 200)
        try:
            ax2[0][lam_tracker].scatter(eta_attractors[rho]['2'][0], plotting_dict_deaths_2[rho]['2'][0], edgecolors = f"{colours[i]}", marker = "s", facecolors= "none", s = 200)
            ax2[0][lam_tracker].scatter(eta_attractors[rho]['2'][-1], plotting_dict_deaths_2[rho]['2'][-1], edgecolors = f"{colours[i]}", marker = "s", facecolors= f"{colours[i]}", s = 200)
        except:
            pass
#         ax[1][lam_tracker].plot(plotting_dict_hyper[rho]["r"], plotting_dict_deaths_2[rho]["r"],  c = f"{colours[i]}", linestyle='dashed')


    ax2[0][lam_tracker].plot(plotting_dict_hyper[rho]["a"], [1 for val in plotting_dict_hyper[rho]["a"]], c = "k")
    ax2[0][lam_tracker].plot(plotting_dict_hyper[rho]["r"], [1 for val in plotting_dict_hyper[rho]["r"]], c = "k")

#     ax[1][lam_tracker].set_xlabel(xlabel, fontsize = 24)
    ax2[0][lam_tracker].set_xlabel("", fontsize = 0)
    ax2[0][lam_tracker].set_ylim([0.8, 3.2])
    ax2[0][lam_tracker].text(0.05,0.95,texts[0][lam_tracker], transform=ax2[0][lam_tracker].transAxes, fontsize = 30)
    ax2[0][lam_tracker].tick_params(axis='both', which='major', labelsize=34)
    if lam_tracker == 0:
        ax2[0][lam_tracker].set_ylabel("Relative Death rate", fontsize = 34)
#     ax[1][lam_tracker].set_title(f"Relative Infected Burden, $\lambda$ = {lam}", fontsize = 20)
    ax2[0][lam_tracker].set_title(fr"{titles[lam_tracker]}($\lambda$ = {lam})", fontsize = 34)
    

    for i, rho in enumerate(rhos):
        ax2[1][lam_tracker].plot(eta_attractors[rho]['1'], plotting_dict_pops[rho]["1"], c = f"{colours[i]}")
        ax2[1][lam_tracker].plot(eta_attractors[rho]['2'], plotting_dict_pops[rho]["2"], c = f"{colours[i]}")
        ax2[1][lam_tracker].scatter(eta_attractors[rho]['1'][0], plotting_dict_pops[rho]['1'][0], edgecolors = f"{colours[i]}", marker = "o", facecolors="none", s = 200)
        ax2[1][lam_tracker].scatter(eta_attractors[rho]['1'][-1], plotting_dict_pops[rho]['1'][-1], edgecolors = f"{colours[i]}", marker = "o", facecolors= f"{colours[i]}", s = 200)
        try:
            ax2[1][lam_tracker].scatter(eta_attractors[rho]['2'][-1], plotting_dict_pops[rho]['2'][-1], edgecolors = f"{colours[i]}", marker = "s", facecolors= f"{colours[i]}", s = 200)
            ax2[1][lam_tracker].scatter(eta_attractors[rho]['2'][0], plotting_dict_pops[rho]['2'][0], edgecolors = f"{colours[i]}", marker = "s", facecolors="none", s = 200)
        except:
            pass
#         ax[2][lam_tracker].plot(plotting_dict_hyper[rho]["r"], plotting_dict_pops[rho]["r"],  c = f"{colours[i]}", linestyle='dashed')

    ax2[1][lam_tracker].plot(plotting_dict_hyper[rho]["a"], [1 for val in plotting_dict_hyper[rho]["a"]], c = "k")
    ax2[1][lam_tracker].plot(plotting_dict_hyper[rho]["r"], [1 for val in plotting_dict_hyper[rho]["r"]], c = "k")


    if lam_tracker == 1:
        ax2[1][lam_tracker].set_xlabel(xlabel, fontsize = 34)
    else:
        ax2[1][lam_tracker].set_xlabel("", fontsize = 0)
        
    ax2[1][lam_tracker].set_ylim([0.3, 1.5])
    ax2[1][lam_tracker].tick_params(axis='both', which='major', labelsize=34)
    ax2[1][lam_tracker].text(0.05,0.95,texts[1][lam_tracker], transform=ax2[1][lam_tracker].transAxes, fontsize = 30)
    
    if lam_tracker == 0:
        ax2[1][lam_tracker].set_ylabel("Relative Population Size", fontsize = 34)
#     ax[2][lam_tracker].set_title(f"Relative Population Size, $\lambda$ = {lam}", fontsize = 20)


fig2.legend(lines2,     # The line objects
   labels=labels,
   title='Probability of cotransmission',  # Title for the legend
   fontsize = 24,
   loc="upper right",
   title_fontsize=26,
   bbox_to_anchor=(0.90,0.85)
   )
# plt.subplots_adjust(right=0.79)
plt.savefig(f"{output_folder}/effects_of_hyperparasitism.png", bbox_inches='tight')
plt.close()

fig1.legend(lines1,     # The line objects
   labels=labels,
   title='Probability of cotransmission',  # Title for the legend
   fontsize = 24,
   loc="upper right",
   title_fontsize=26,
   bbox_to_anchor=(0.90,0.85)
   )
# plt.subplots_adjust(right=0.79)
plt.savefig(f"{output_folder}/evolved_virulences.png", bbox_inches='tight')
plt.close()
del sol