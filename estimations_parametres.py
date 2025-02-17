#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:59:35 2024

@author: alexandreperrin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy import stats, sparse
import pandas as pd
from scipy.optimize import minimize
from matplotlib.colors import LinearSegmentedColormap
import ot
import cma

# data_si_mops_acet = pd.read_excel("1-s2.0-S0960982219304919-mmc3.xlsx",sheet_name = 'MG1655 M9 acetate', header=None)
# data_si_mops_acet = np.array(data_si_mops_acet)
# data_si_mops_glu1 = pd.read_excel("1-s2.0-S0960982219304919-mmc3.xlsx",sheet_name = 'MG1655 MOPS glucose', header=None)
# data_si_mops_glu1 = np.array(data_si_mops_glu1)
# data_si_mops_gly11aa = pd.read_excel("1-s2.0-S0960982219304919-mmc3.xlsx",sheet_name = 'MG1655 MOPS glycerol 11aa', header=None)
# data_si_mops_gly11aa = np.array(data_si_mops_gly11aa)
# data_si_mops_arg = pd.read_excel("1-s2.0-S0960982219304919-mmc3.xlsx",sheet_name = 'NCM3722 MOPS arginine', header=None)
# data_si_mops_arg = np.array(data_si_mops_arg)
# data_si_mops_glu12aa = pd.read_excel("1-s2.0-S0960982219304919-mmc3.xlsx",sheet_name = 'NCM3722 MOPS glucose 12aa', header=None)
# data_si_mops_glu12aa = np.array(data_si_mops_glu12aa)
# data_si_mops_glu2 = pd.read_excel("1-s2.0-S0960982219304919-mmc3.xlsx",sheet_name = 'NCM3722 MOPS glucose', header=None)
# data_si_mops_glu2 = np.array(data_si_mops_glu2)
# data_si_mops_glu1 = np.delete(data_si_mops_glu1, 308, axis = 0)

# data_Tiruvadi_JM85_GluCas = pd.read_csv("Fig1-2_Data/JM85(ssb-mYpet)_GluCas.txt", sep = '\t', header = None)
# data_Tiruvadi_JM85_GluCas = np.array(data_Tiruvadi_JM85_GluCas)
# data_Tiruvadi_JM85_Gly = pd.read_csv("Fig1-2_Data/JM85(ssb-mYpet)_Gly.txt", sep = '\t', header = None)
# data_Tiruvadi_JM85_Gly = np.array(data_Tiruvadi_JM85_Gly)
# data_Tiruvadi_STK13_Ace = pd.read_csv("Fig1-2_Data/STK13(DnaN-FtsN)_Ace.txt", sep = '\t', header = None)
# data_Tiruvadi_STK13_Ace = np.array(data_Tiruvadi_STK13_Ace)
# data_Tiruvadi_STK13_AlaTrE = pd.read_csv("Fig1-2_Data/STK13(DnaN-FtsN)_AlaTrE.txt", sep = '\t', header = None)
# data_Tiruvadi_STK13_AlaTrE = np.array(data_Tiruvadi_STK13_AlaTrE)
# data_Tiruvadi_STK13_Glu = pd.read_csv("Fig1-2_Data/STK13(DnaN-FtsN)_Glu.txt", sep = '\t', header = None)
# data_Tiruvadi_STK13_Glu = np.array(data_Tiruvadi_STK13_Glu)
# data_Tiruvadi_STK13_GluCas= pd.read_csv("Fig1-2_Data/STK13(DnaN-FtsN)_GluCas.txt", sep = '\t', header = None)
# data_Tiruvadi_STK13_GluCas = np.array(data_Tiruvadi_STK13_GluCas)
# data_Tiruvadi_STK13_Gly= pd.read_csv("Fig1-2_Data/STK13(DnaN-FtsN)_Gly.txt", sep = '\t', header = None)
# data_Tiruvadi_STK13_Gly = np.array(data_Tiruvadi_STK13_Gly)
# data_Tiruvadi_STK13_GlyCas= pd.read_csv("Fig1-2_Data/STK13(DnaN-FtsN)_GlyCas.txt", sep = '\t', header = None)
# data_Tiruvadi_STK13_GlyCas = np.array(data_Tiruvadi_STK13_GlyCas)
# data_Tiruvadi_STK13_GlyTrE= pd.read_csv("Fig1-2_Data/STK13(DnaN-FtsN)_GlyTrE.txt", sep = '\t', header = None)
# data_Tiruvadi_STK13_GlyTrE = np.array(data_Tiruvadi_STK13_GlyTrE)
# data_Tiruvadi_STK13_Man= pd.read_csv("Fig1-2_Data/STK13(DnaN-FtsN)_Man.txt", sep = '\t', header = None)
# data_Tiruvadi_STK13_Man = np.array(data_Tiruvadi_STK13_Man)

data_witz = pd.read_csv('Fig1_2_3.csv', header = None) 
[Lb, dL, Li, Li_prev, dLambda_if, lambda_inv, Tid, dLambda_id, Lb_mother, condition] = range(0, 10)



    





    




def extract_donnees(data_or_cdt, author):
    if author == 'si':
        data = data_or_cdt
        v = data[1:, data[0] == 'newborn size (micron)'].astype(float)
        v_tild = data[1:, data[0] == 'initiation size per ori (micron)'].astype(float)/2.
        Delta_app = data[1:, data[0] == 'added size (micron)'].astype(float)
        R_app = data[1:, data[0] == 'tau_cyc (minute)'].astype(float)/60.
        V = (Delta_app + v)/2.
        lamb = np.log(V/v_tild)/R_app
        return [v, v_tild, V, Delta_app, R_app, lamb]
    
    elif author == 'tiru':
        data = data_or_cdt
        v = data[1:, data[0] == 'Lb'].astype(float)
        v_tild = data[1:, data[0] == 'Lri'].astype(float)/2.
        Vd = data[1:, data[0] == 'Ld'].astype(float)
        Tri = data[1:, data[0] == 'Tri'].astype(float)
        Td = data[1:, data[0] == 'Td'].astype(float)
        R_app = (Td - Tri)/60.
        Delta_app = Vd - v
        V = Vd/2.
        lamb = np.log(V/v_tild)/R_app
        return [v, v_tild, V, Delta_app, R_app, lamb]
    
    elif author == 'witz':
        data = data_witz[data_witz[condition] == data_or_cdt]
        data = data.dropna(subset = (Lb, dL, Li, Tid))
        data = data.to_numpy()
        data = data[:,:-1].astype(float)
        Delta_app = data[:, dL]
        R_app = data[:, Tid]/60.
        v = data[:, Lb]
        V = (v + Delta_app)/2.
        v_tild = data[:, Li]/2.
        lamb = np.log(V/v_tild)/R_app
        return [v, v_tild, V, Delta_app, R_app, lamb]
    
def extract_donnees2(data_or_cdt, author):
    if author == 'si':
        data = data_or_cdt
        v = data[1:, data[0] == 'newborn size (micron)'].astype(float)
        v_tild = data[1:, data[0] == 'initiation size per ori (micron)'].astype(float)/2.
        V = data[1:, data[0] == 'division size (micron)'].astype(float)
        C = data[1:, data[0] == 'C period (minute)'].astype(float)/60.
        tau = data[1:, data[0] == 'generation time (minute)'].astype(float)/60.
        lamb = np.log(V/v)/tau
        return [v, V, v_tild, C, lamb]
    
    elif author == 'tiru':
        data = data_or_cdt
        v = data[1:, data[0] == 'Lb'].astype(float)
        v_tild = data[1:, data[0] == 'Lri'].astype(float)/2.
        V = data[1:, data[0] == 'Ld'].astype(float)
        Tri = data[1:, data[0] == 'Tri'].astype(float)
        Trt = data[1:, data[0] == 'Trt'].astype(float)
        tau = data[1:, data[0] == 'Td'].astype(float)
        C = Trt - Tri
        lamb = np.log(V/v)/tau
        return [v, V, v_tild, C, lamb]
    
def extract_donnees3(data_or_cdt, author):
    if author == 'si':
        I = adder_init(data_or_cdt)
        Delta_i = I[1] - I[0]/2.
        data = data_or_cdt
        v_tild = data[1:, data[0] == 'initiation size per ori (micron)'].astype(float)/2.
        V = data[1:, data[0] == 'division size (micron)'].astype(float)/2.
        return [V - v_tild, Delta_i]
    
    elif author == 'witz':
        data = data_witz[data_witz[condition] == data_or_cdt]
        data = data.dropna(subset = (dLambda_id, dLambda_if))
        data = data.to_numpy()
        data = data[:,:-1].astype(float)
        Delta_i = data[:, dLambda_if]
        Delta_id = data[:, dLambda_id]/2.
        return [Delta_id, Delta_i]
    

    
def extract_donnees32(data_or_cdt, author):
    if author == 'si':
        I = adder_init(data_or_cdt)
        Delta_i = I[1] - I[0]/2.
        data = data_or_cdt
        v = data[1:, data[0] == 'newborn size (micron)'].astype(float)
        V = data[1:, data[0] == 'division size (micron)'].astype(float)
        Delta_d = V-v
        return [Delta_d.T[0], Delta_i]
    
    elif author == 'witz':
        data = data_witz[data_witz[condition] == data_or_cdt]
        data = data.dropna(subset = (dLambda_id, dLambda_if))
        data = data.to_numpy()
        data = data[:,:-1].astype(float)
        Delta_i = data[:, dLambda_if]
        Delta_d = data[:, dL]
        return [Delta_d, Delta_i]
    
def figure_overlap_Dd_Di(data_or_cdt, author):
    [Delta_d, Delta_i] = extract_donnees32(data_or_cdt, author)
    binn = np.linspace(0, 5, 100)
    plt.hist(Delta_d, binn, color = 'k', alpha = 0.6, density = True, label = '$\Delta_d$')
    plt.hist(Delta_i, binn, color = 'r', alpha = 0.6, density = True, label = '$\Delta_i$')
    plt.legend()
    plt.ylabel("$Density$")
    plt.xlabel("$Size$ ($\mu m$)")
    

def extract_donnees4(data_or_cdt, author):
    if author == 'si':
        data = select_data_si_seq_vi(data_or_cdt)
        vb_1 = data[1:, data[0] == 'previous newborn size (micron)'].astype(float)
        vb_2 =vb_1 = data[1:, data[0] == 'newborn size (micron)'].astype(float)
        vi_1 = data[1:, data[0] == 'previous initiation size per ori (micron)'].astype(float)/2.
        vi_2 = data[1:, data[0] == 'initiation size per ori (micron)'].astype(float)/2.
        CD_1 = data[1:, data[0] == 'previous tau_cyc (minute)'].astype(float)
        lamb_1 = np.log(vb_2/vi_1)/CD_1
        return [vb_1, vi_1, vi_2, lamb_1]
    
    elif author == 'witz':
        data = data_witz[data_witz[condition] == data_or_cdt]
        data = data.dropna(subset = (Li, Li_prev, Lb_mother, dLambda_if))
        data = data.to_numpy()
        data = data[:,:-1].astype(float)
        vi_1 = data[:, Li_prev]
        vi_2 = data[:, Li]/2.
        vb_1 = data[:, Lb_mother]
        return [vb_1, vi_1, vi_2]

        
    
def MnteCarlo_proba_DIAM(data_or_cdt, author):
    [Delta_d, Delta_i] = extract_donnees32(data_or_cdt, author)
    n = 1000000
    Vb = np.random.choice(Delta_d, size = n)
    Vi = np.random.choice(Delta_i, size = n)
    for loop in range(100):
        Vb = (Vb + np.random.choice(Delta_d, size = n))/2.
        Vi = (Vi + np.random.choice(Delta_i, size = n))/2.
    Vb = np.reshape(Vb, (1, n))
    Vi = np.reshape(Vi, (1, n))
    Vb = np.concatenate((Vb, np.ones((1, n))), axis = 0)
    Vi = np.concatenate((Vi, np.zeros((1, n))), axis = 0)
    Total = np.concatenate((Vb, Vi), axis = 1)
    i = np.argsort(Total[0])
    Total = Total[:,i]
    s = np.sum(np.arange(1, 2*n + 1)*Total[1])
    plt.hist(Vb[0], [0.01*k for k in range(400)], color = 'b', density = True, alpha = 0.5)
    plt.hist(Vi[0], [0.01*k for k in range(400)], color = 'k', density = True, alpha = 0.5)
    return 1 + n*(n+1)/(2 * n**2) - s/n**2 # Wilcoxon strategy


def estim_proba_Tbdecroiss_empir(data_or_cdt, author):
    [Delta_id, Delta_i] = extract_donnees3(data_or_cdt, author)
    Delta_id = Delta_id.flatten().astype(float)
    Delta_i = Delta_i.flatten().astype(float)
    cpt = 0
    N_id = np.size(Delta_id)
    N_i = np.size(Delta_i)
    for d_id1 in Delta_id:
        for d_id2 in Delta_id:
            for d_i in Delta_i:
                if 2*d_id1 + d_i < d_id2:
                    cpt += 1
    return [cpt/(N_id*N_id*N_i), N_id, N_i]


def estim_param_AND(data_or_cdt, author):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    def f(X):
        [a_d, s_d, a_r, s_r] = X
        return - np.sum(np.log(2*stats.gamma.pdf(Delta_app, a = a_d, scale = s_d)*stats.gamma.cdf(R_app, a = a_r, scale = s_r) + stats.gamma.cdf(Delta_app, a = a_d, scale = s_d)*stats.gamma.pdf(R_app, a = a_r, scale = s_r)/(lamb*V)))
    moy_d = np.mean(Delta_app)
    var_d = np.var(Delta_app)
    moy_r = np.mean(R_app)
    var_r = np.var(R_app)
    X0 = [moy_d**2/var_d, var_d/moy_d, moy_r**2/var_r, var_r/moy_r]

    results = minimize(f, X0, method = 'Nelder-Mead')
    BIC = results.fun + 2*np.log(np.size(v))
    return [results.x, results.success, results.message, BIC]

def estim_param_revolution(data_or_cdt, author):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    def f(X):
        [a_d, s_d, a_r, s_r] = X
        return - np.sum(np.log(2*stats.gamma.pdf(Delta_app, a = a_d, scale = s_d)*stats.gamma.cdf(V - v_tild, a = a_r, scale = s_r) + stats.gamma.cdf(Delta_app, a = a_d, scale = s_d)*stats.gamma.pdf(V - v_tild, a = a_r, scale = s_r)))
    moy_d = np.mean(Delta_app)
    var_d = np.var(Delta_app)
    moy_id = np.mean(V - v_tild)
    var_id = np.var(V - v_tild)
    X0 = [moy_d**2/var_d, var_d/moy_d, moy_id**2/var_id, var_id/moy_id]

    results = minimize(f, X0, method = 'Nelder-Mead')
    BIC = results.fun + 2*np.log(np.size(v))
    return [results.x, results.success, results.message, BIC]

def estim_param_revolution_test_densites(data_or_cdt, author, law_Delta_d, law_Delta_id):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    def f(X):
        [a_d, s_d, a_r, s_r] = X
        return - np.sum(np.log(2*law_Delta_d.pdf(Delta_app, a_d, scale = s_d)*law_Delta_id.cdf(V - v_tild, a_r, scale = s_r) + law_Delta_d.cdf(Delta_app, a_d, scale = s_d)*law_Delta_id.pdf(V - v_tild, a_r, scale = s_r)))
    moy_d = np.mean(Delta_app)
    var_d = np.var(Delta_app)
    moy_id = np.mean(V - v_tild)
    var_id = np.var(V - v_tild)
    if law_Delta_d == stats.gamma:
        a_d = moy_d**2/var_d
        scale_d = var_d/moy_d
        if law_Delta_id == stats.gamma:
            a_id = moy_id**2/var_id
            scale_id = var_id/moy_id
            X0 = [a_d, scale_d, a_id, scale_id]
        elif law_Delta_id == stats.lognorm:
            s_id = np.sqrt(np.log(var_id/moy_id**2 + 1))
            scale_id = moy_id
            X0 = [a_d, scale_d, s_id, scale_id]
    elif law_Delta_d == stats.lognorm:
        a_d = np.sqrt(np.log(var_d/moy_d**2 + 1))
        scale_d = moy_d
        if law_Delta_id == stats.gamma:
            a_id = moy_id**2/var_id
            scale_id = var_id/moy_id
            X0 = [a_d, scale_d, a_id, scale_id]
        elif law_Delta_id == stats.lognorm:
            s_id = np.sqrt(np.log(var_id/moy_id**2 + 1))
            scale_id = moy_id
            X0 = [a_d, scale_d, s_id, scale_id]
    results = cma.fmin(f, X0, sigma0 = 1, options = {'verbose': -9})
    BIC = results[1] + 2*np.log(np.size(v))
    return [results[0], BIC]


def likelihhod_revolution(law_Delta_d, law_Delta_id):
    DATA = [data_Tiruvadi_STK13_Ace, data_Tiruvadi_STK13_AlaTrE, data_si_mops_acet, data_Tiruvadi_JM85_Gly, data_Tiruvadi_STK13_Gly, data_Tiruvadi_STK13_GlyTrE, data_Tiruvadi_STK13_Glu, data_si_mops_arg, 'glycerol', data_Tiruvadi_STK13_GlyCas, data_si_mops_gly11aa, 'glucose', data_si_mops_glu1, 'glucose8a',  data_si_mops_glu2, data_si_mops_glu12aa]
    AUTHOR = ['tiru', 'tiru', 'si', 'tiru', 'tiru', 'tiru', 'tiru', 'si', 'witz', 'tiru', 'si', 'witz', 'si', 'witz', 'si', 'si']
    for i in range(len(DATA)):
        print(estim_param_revolution_test_densites(DATA[i],AUTHOR[i], law_Delta_d, law_Delta_id))
    

def figure_simu_CPM(data_or_cdt, author, nb_sim, alpha):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    v_tild = v_tild.flatten()
    V = V.flatten()
    v = v.flatten()
    lamb = lamb.flatten()
    mean_V = np.mean(V)
    sd_V = np.std(V)
    N = np.size(v)
    maxX = mean_V + 3*sd_V
    minX = mean_V - 2.5*sd_V
    
    X = np.linspace(minX, maxX, 200)
    [a_d, s_d, a_r, s_r]  = estim_param_AND(data_or_cdt, author)[0]
    hatF_Vbsim = np.zeros((nb_sim, 200))
    for n in range(nb_sim):
        R = stats.gamma.rvs(size = N, a = a_r, scale = s_r)
        Delta_d = stats.gamma.rvs(size = N, a = a_d, scale = s_d)
        Rep = v_tild*np.exp(lamb*R)
        Div = (v + Delta_d)/2.
        cond_rep = np.where(Rep > Div, 1, 0)
        Vb_sim = Rep*cond_rep + Div*(1 - cond_rep)
        hatF_Vbsim[n, :] = funct_rep1(Vb_sim, X)
    hatF_Vbsim = np.sort(hatF_Vbsim, axis = 0)
    bornesup = hatF_Vbsim[int((1-alpha/2.)*nb_sim), :]
    borneinf = hatF_Vbsim[int(alpha*nb_sim/2.), :]
    Y = funct_rep1(V, X)
    plt.fill_between(X, borneinf, bornesup, color = 'slateblue', alpha = 0.5)
    plt.plot(X, Y, 'k')
    
def figure_simu_CPM2(data_or_cdt, author, nb_sim, alpha):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    v_tild = v_tild.flatten()
    V = V.flatten()
    v = v.flatten()
    lamb = lamb.flatten()
    E_V = np.mean(V)
    E_v_tild = np.mean(v_tild)
    N = np.size(v)
    Xvb = np.linspace(0, 2*E_V, 200)
    Xvi = np.linspace(0, 2*E_v_tild, 200)
    [a_d, s_d, a_r, s_r]  = estim_param_AND(data_or_cdt, author)[0]
    hatF_Vbsim = np.zeros((nb_sim, 200, 200))
    for n in range(nb_sim):
        R = stats.gamma.rvs(size = N, a = a_r, scale = s_r)
        Delta_d = stats.gamma.rvs(size = N, a = a_d, scale = s_d)
        Rep = v_tild*np.exp(lamb*R)
        Div = (v + Delta_d)/2.
        cond_rep = np.where(Rep > Div, 1, 0)
        Vb_sim = Rep*cond_rep + Div*(1 - cond_rep)
        hatF_Vbsim[n, :, :] = funct_rep(Vb_sim, v_tild, Xvb, Xvi)
    hatF_Vbsim = np.sort(hatF_Vbsim, axis = 0)
    meanF = np.mean(hatF_Vbsim , axis = 0)
    bornesup = hatF_Vbsim[int((1-alpha/2.)*nb_sim), :, :] - meanF
    borneinf = hatF_Vbsim[int(alpha*nb_sim/2.), :, :] - meanF
    F_exp = funct_rep(V, v_tild, Xvb, Xvi)
    Dif = F_exp - meanF
    cond_over = np.where(F_exp > bornesup, 1, 0)
    cond_under = np.where(F_exp < borneinf, 1, 0)
    out = Dif*(cond_over + cond_under)
    X, Y = np.meshgrid(Xvi, Xvb)
    plt.pcolormesh(X, Y, out, cmap=custom_cmap, vmin=-0.05, vmax=0.05)

def figure_simu_DIAM(data_or_cdt, author, nb_sim, alpha):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    V = V.flatten()
    v = v.flatten()
    mean_V = np.mean(V)
    sd_V = np.std(V)
    N = np.size(v)
    maxX = mean_V + 3*sd_V
    minX = mean_V - 2.5*sd_V
    X = np.linspace(minX, maxX, 200)
    [a, c, scale, BIC]  = estim_param_DIAM(data_or_cdt, author)
    hatF_Vbsim = np.zeros((nb_sim, 200))
    for n in range(nb_sim):
        Delta_d = stats.gengamma.rvs(size = N, a = a, c = c, scale = scale)
        Vb_sim = (v + Delta_d)/2.
        hatF_Vbsim[n, :] = funct_rep1(Vb_sim, X)
    hatF_Vbsim = np.sort(hatF_Vbsim, axis = 0)
    bornesup = hatF_Vbsim[int((1-alpha/2.)*nb_sim), :]
    borneinf = hatF_Vbsim[int(alpha*nb_sim/2.), :]
    Y = funct_rep1(V, X)
    plt.fill_between(X, borneinf, bornesup, color = 'mediumaquamarine', alpha = 0.5)
    plt.plot(X, Y, 'k')

def figure_simu_RDAM(data_or_cdt, author, nb_sim, alpha):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    V = V.flatten()
    v = v.flatten()
    mean_V = np.mean(V)
    sd_V = np.std(V)
    N = np.size(v)
    maxX = mean_V + 3*sd_V
    minX = mean_V - 2.5*sd_V
    X = np.linspace(minX, maxX, 200)
    [a, c, scale, BIC]  = estim_param_RDAM(data_or_cdt, author)
    hatF_Vbsim = np.zeros((nb_sim, 200))
    for n in range(nb_sim):
        Delta_id = stats.gengamma.rvs(size = N, a = a, c = c, scale = scale)
        Vb_sim = v_tild + Delta_id
        hatF_Vbsim[n, :] = funct_rep1(Vb_sim, X)
    hatF_Vbsim = np.sort(hatF_Vbsim, axis = 0)
    bornesup = hatF_Vbsim[int((1-alpha/2.)*nb_sim), :]
    borneinf = hatF_Vbsim[int(alpha*nb_sim/2.), :]
    Y = funct_rep1(V, X)
    plt.fill_between(X, borneinf, bornesup, color = 'darkcyan', alpha = 0.5)
    plt.plot(X, Y, 'k')
    


def figure_sim1(data, author, nbsim, alpha):
    plt.subplot(131)
    figure_simu_DIAM(data, author, nbsim, alpha)
    plt.subplot(132)
    figure_simu_CPM(data, author, nbsim, alpha)
    plt.subplot(133)
    figure_simu_RDAM(data, author, nbsim, alpha)
    

def figure_simu_DIAM2(data_or_cdt, author, nb_sim, alpha):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    V = V.flatten()
    v = v.flatten()
    v_tild = v_tild.flatten()
    E_V = np.mean(V)
    E_v_tild = np.mean(v_tild)
    N = np.size(v)
    Xvb = np.linspace(0, 2*E_V, 200)
    Xvi = np.linspace(0, 2*E_v_tild, 200)
    [a, c, scale, BIC]  = estim_param_DIAM(data_or_cdt, author)
    hatF_Vbsim = np.zeros((nb_sim, 200, 200))
    for n in range(nb_sim):
        Delta_d = stats.gengamma.rvs(size = N, a = a, c = c, scale = scale)
        Vb_sim = (v + Delta_d)/2.
        hatF_Vbsim[n, :, :] = funct_rep(Vb_sim, v_tild, Xvb, Xvi)
    hatF_Vbsim = np.sort(hatF_Vbsim, axis = 0)
    meanF = np.mean(hatF_Vbsim , axis = 0)
    bornesup = hatF_Vbsim[int((1-alpha/2.)*nb_sim), :, :] - meanF
    borneinf = hatF_Vbsim[int(alpha*nb_sim/2.), :, :] - meanF
    F_exp = funct_rep(V, v_tild, Xvb, Xvi)
    Dif = F_exp - meanF
    cond_over = np.where(F_exp > bornesup, 1, 0)
    cond_under = np.where(F_exp < borneinf, 1, 0)
    out = Dif*(cond_over + cond_under)
    X, Y = np.meshgrid(Xvi, Xvb)
    plt.pcolormesh(X, Y, out, cmap=custom_cmap, vmin=-0.05, vmax=0.05)

def figure_simu_RDAM2(data_or_cdt, author, nb_sim, alpha):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    V = V.flatten()
    v_tild = v_tild.flatten()
    E_V = np.mean(V)
    E_v_tild = np.mean(v_tild)
    N = np.size(v)
    Xvb = np.linspace(0, 2*E_V, 200)
    Xvi = np.linspace(0, 2*E_v_tild, 200)
    [a, c, scale, BIC]  = estim_param_RDAM(data_or_cdt, author)
    hatF_Vbsim = np.zeros((nb_sim, 200, 200))
    for n in range(nb_sim):
        Delta_id = stats.gengamma.rvs(size = N, a = a, c = c, scale = scale)
        Vb_sim = v_tild + Delta_id
        hatF_Vbsim[n, :, :] = funct_rep(Vb_sim, v_tild, Xvb, Xvi)
    hatF_Vbsim = np.sort(hatF_Vbsim, axis = 0)
    meanF = np.mean(hatF_Vbsim , axis = 0)
    bornesup = hatF_Vbsim[int((1-alpha/2.)*nb_sim), :, :] - meanF
    borneinf = hatF_Vbsim[int(alpha*nb_sim/2.), :, :] - meanF
    F_exp = funct_rep(V, v_tild, Xvb, Xvi)
    Dif = F_exp - meanF
    cond_over = np.where(F_exp > bornesup, 1, 0)
    cond_under = np.where(F_exp < borneinf, 1, 0)
    out = Dif*(cond_over + cond_under)

    X, Y = np.meshgrid(Xvi, Xvb)
    plt.pcolormesh(X, Y, out, cmap=custom_cmap, vmin=-0.05, vmax=0.05)
    plt.colorbar()
    
custom_cmap = LinearSegmentedColormap.from_list('PurpleWhiteCyan', colors = ('darkcyan', 'w', 'slateblue'), N=256)    


def figure_simu_RDAM20(data_or_cdt, author, nb_sim, alpha, index):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    N = np.size(v)
    V = np.reshape(V, (1,N))
    v_tild = np.reshape(v_tild, (1,N))
    v = np.reshape(v, (1,N))
    lamb = np.reshape(lamb, (1,N))
    A = np.concatenate((V, v_tild, v, lamb), axis = 0)
    sorted_indices = np.argsort(A[index, :])
    A = A[:, sorted_indices]
    vb = A[0,:]
    vi = A[1,:]
    Nsub = int(N/4.)
    
    mean_vb = np.mean(vb)
    sd_vb = np.std(vb)

    maxXvb = mean_vb + 3*sd_vb
    minXvb = mean_vb - 2.5*sd_vb

    Xvb = np.linspace(minXvb, maxXvb, 200)
    [a, c, scale, BIC]  = estim_param_RDAM(data_or_cdt, author)
    hatF_Vbsim = np.zeros((nb_sim, 4, 200))
    for n in range(nb_sim):
        Delta_id = stats.gengamma.rvs(size = N, a = a, c = c, scale = scale)
        Vb_sim = vi + Delta_id
        hatF_Vbsim[n, :, :] = funct_rep2(Vb_sim, Nsub, Xvb)
    hatF_Vbsim = np.sort(hatF_Vbsim, axis = 0)
    bornesup = hatF_Vbsim[int((1-alpha/2.)*nb_sim), :, :]
    borneinf = hatF_Vbsim[int(alpha*nb_sim/2.), :, :]
    F_exp = funct_rep2(vb, Nsub, Xvb)
    plt.subplot(4,4,9)
    plt.fill_between(Xvb, borneinf[0,:], bornesup[0,:], color = 'darkcyan', alpha = 0.5)
    plt.plot(Xvb, F_exp[0,:], 'k')
    plt.subplot(4,4,10)
    plt.fill_between(Xvb, borneinf[1,:], bornesup[1,:], color = 'darkcyan', alpha = 0.5)
    plt.plot(Xvb, F_exp[1,:], 'k')
    plt.subplot(4,4,11)
    plt.fill_between(Xvb, borneinf[2,:], bornesup[2,:], color = 'darkcyan', alpha = 0.5)
    plt.plot(Xvb, F_exp[2,:], 'k')
    plt.subplot(4,4,12)
    plt.fill_between(Xvb, borneinf[3,:], bornesup[3,:], color = 'darkcyan', alpha = 0.5)
    plt.plot(Xvb, F_exp[3,:], 'k')  

def figure_simu_CPM20(data_or_cdt, author, nb_sim, alpha, index):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    N = np.size(v)
    V = np.reshape(V, (1,N))
    v_tild = np.reshape(v_tild, (1,N))
    v = np.reshape(v, (1,N))
    lamb = np.reshape(lamb, (1,N))
    A = np.concatenate((V, v_tild, v, lamb), axis = 0)
    sorted_indices = np.argsort(A[index, :])
    A = A[:, sorted_indices]
    vb2 = A[0,:]
    vb1 = A[2,:]
    vi = A[1,:]
    lamb = A[3,:]
    Nsub = int(N/4.)
    
    mean_vb = np.mean(vb2)
    sd_vb = np.std(vb2)

    maxXvb = mean_vb + 3*sd_vb
    minXvb = mean_vb - 2.5*sd_vb

    Xvb = np.linspace(minXvb, maxXvb, 200)
    [a_d, s_d, a_r, s_r]  = estim_param_AND(data_or_cdt, author)[0]
    hatF_Vbsim = np.zeros((nb_sim, 4, 200))
    for n in range(nb_sim):
        R = stats.gamma.rvs(size = N, a = a_r, scale = s_r)
        Delta_d = stats.gamma.rvs(size = N, a = a_d, scale = s_d)
        Rep = vi*np.exp(lamb*R)
        Div = (vb1 + Delta_d)/2.
        cond_rep = np.where(Rep > Div, 1, 0)
        Vb_sim = Rep*cond_rep + Div*(1 - cond_rep)
        hatF_Vbsim[n, :, :] = funct_rep2(Vb_sim, Nsub, Xvb)
    hatF_Vbsim = np.sort(hatF_Vbsim, axis = 0)
    bornesup = hatF_Vbsim[int((1-alpha/2.)*nb_sim), :, :]
    borneinf = hatF_Vbsim[int(alpha*nb_sim/2.), :, :]
    F_exp = funct_rep2(vb2, Nsub, Xvb)
    plt.subplot(4,4,5)
    plt.fill_between(Xvb, borneinf[0,:], bornesup[0,:], color = 'slateblue', alpha = 0.5)
    plt.plot(Xvb, F_exp[0,:], 'k')
    plt.subplot(4,4,6)
    plt.fill_between(Xvb, borneinf[1,:], bornesup[1,:], color = 'slateblue', alpha = 0.5)
    plt.plot(Xvb, F_exp[1,:], 'k')
    plt.subplot(4,4,7)
    plt.fill_between(Xvb, borneinf[2,:], bornesup[2,:], color = 'slateblue', alpha = 0.5)
    plt.plot(Xvb, F_exp[2,:], 'k')
    plt.subplot(4,4,8)
    plt.fill_between(Xvb, borneinf[3,:], bornesup[3,:], color = 'slateblue', alpha = 0.5)
    plt.plot(Xvb, F_exp[3,:], 'k')
    
    
def figure_simu_CAM20(data_or_cdt, author, nb_sim, alpha, index):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    N = np.size(v)
    V = np.reshape(V, (1,N))
    v_tild = np.reshape(v_tild, (1,N))
    v = np.reshape(v, (1,N))
    lamb = np.reshape(lamb, (1,N))
    A = np.concatenate((V, v_tild, v, lamb), axis = 0)
    sorted_indices = np.argsort(A[index, :])
    A = A[:, sorted_indices]
    vb2 = A[0,:]
    vb1 = A[2,:]
    vi = A[1,:]
    lamb = A[3,:]
    Nsub = int(N/4.)
    
    mean_vb = np.mean(vb2)
    sd_vb = np.std(vb2)

    maxXvb = mean_vb + 3*sd_vb
    minXvb = mean_vb - 2.5*sd_vb

    Xvb = np.linspace(minXvb, maxXvb, 200)
    [a_d, s_d, a_r, s_r]  = estim_param_revolution(data_or_cdt, author)[0]
    hatF_Vbsim = np.zeros((nb_sim, 4, 200))
    for n in range(nb_sim):
        Delta_id = stats.gamma.rvs(size = N, a = a_r, scale = s_r)
        Delta_d = stats.gamma.rvs(size = N, a = a_d, scale = s_d)
        Rep = vi + Delta_id
        Div = (vb1 + Delta_d)/2.
        cond_rep = np.where(Rep > Div, 1, 0)
        Vb_sim = Rep*cond_rep + Div*(1 - cond_rep)
        hatF_Vbsim[n, :, :] = funct_rep2(Vb_sim, Nsub, Xvb)
    hatF_Vbsim = np.sort(hatF_Vbsim, axis = 0)
    bornesup = hatF_Vbsim[int((1-alpha/2.)*nb_sim), :, :]
    borneinf = hatF_Vbsim[int(alpha*nb_sim/2.), :, :]
    F_exp = funct_rep2(vb2, Nsub, Xvb)
    plt.subplot(4,4,13)
    plt.fill_between(Xvb, borneinf[0,:], bornesup[0,:], color = 'darkred', alpha = 0.5)
    plt.plot(Xvb, F_exp[0,:], 'k')
    plt.subplot(4,4,14)
    plt.fill_between(Xvb, borneinf[1,:], bornesup[1,:], color = 'darkred', alpha = 0.5)
    plt.plot(Xvb, F_exp[1,:], 'k')
    plt.subplot(4,4,15)
    plt.fill_between(Xvb, borneinf[2,:], bornesup[2,:], color = 'darkred', alpha = 0.5)
    plt.plot(Xvb, F_exp[2,:], 'k')
    plt.subplot(4,4,16)
    plt.fill_between(Xvb, borneinf[3,:], bornesup[3,:], color = 'darkred', alpha = 0.5)
    plt.plot(Xvb, F_exp[3,:], 'k')

def figure_simu_DIAM20(data_or_cdt, author, nb_sim, alpha, index):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    N = np.size(v)
    V = np.reshape(V, (1,N))
    v_tild = np.reshape(v_tild, (1,N))
    v = np.reshape(v, (1,N))
    lamb = np.reshape(lamb, (1,N))
    A = np.concatenate((V, v_tild, v, lamb), axis = 0)
    sorted_indices = np.argsort(A[index, :])
    A = A[:, sorted_indices]
    vb2 = A[0,:]
    vb1 = A[2,:]
    Nsub = int(N/4.)

    
    mean_vb = np.mean(vb2)
    sd_vb = np.std(vb2)

    maxXvb = mean_vb + 3*sd_vb
    minXvb = mean_vb - 2.5*sd_vb

    Xvb = np.linspace(minXvb, maxXvb, 200)
    [a, c, scale, BIC]  = estim_param_DIAM(data_or_cdt, author)
    hatF_Vbsim = np.zeros((nb_sim, 4, 200))
    for n in range(nb_sim):
        Delta_d = stats.gengamma.rvs(size = N, a = a, c = c, scale = scale)
        Vb_sim = (vb1 + Delta_d)/2.
        hatF_Vbsim[n, :, :] = funct_rep2(Vb_sim, Nsub, Xvb)
    hatF_Vbsim = np.sort(hatF_Vbsim, axis = 0)
    bornesup = hatF_Vbsim[int((1-alpha/2.)*nb_sim), :, :]
    borneinf = hatF_Vbsim[int(alpha*nb_sim/2.), :, :]
    F_exp = funct_rep2(vb2, Nsub, Xvb)
    plt.subplot(4,4,1)
    plt.fill_between(Xvb, borneinf[0,:], bornesup[0,:], color = 'palevioletred', alpha = 0.5)
    plt.plot(Xvb, F_exp[0,:], 'k')
    plt.subplot(4,4,2)
    plt.fill_between(Xvb, borneinf[1,:], bornesup[1,:], color = 'palevioletred', alpha = 0.5)
    plt.plot(Xvb, F_exp[1,:], 'k')
    plt.subplot(4,4,3)
    plt.fill_between(Xvb, borneinf[2,:], bornesup[2,:], color = 'palevioletred', alpha = 0.5)
    plt.plot(Xvb, F_exp[2,:], 'k')
    plt.subplot(4,4,4)
    plt.fill_between(Xvb, borneinf[3,:], bornesup[3,:], color = 'palevioletred', alpha = 0.5)
    plt.plot(Xvb, F_exp[3,:], 'k')
    
def figure_simu_initadder(data_or_cdt, author, nb_sim, alpha, index):
        [vb_1, vi_1, vi_2, lamb_1] = extract_donnees4(data_or_cdt, author)
        c = ['#ea638c', '#735dbd', '#8798D7']
        N = np.size(vb_1)

        vb_1 = np.reshape(vb_1, (1, N))
        vi_1 = np.reshape(vi_1, (1, N))
        vi_2 = np.reshape(vi_2, (1, N))
        lamb_1 = np.reshape(lamb_1, (1, N))
        A = np.concatenate((vb_1, vi_1, lamb_1, vi_2), axis = 0)
        sorted_indices = np.argsort(A[index, :])
        A = A[:, sorted_indices]
        vi1 = A[1,:]
        vi2 = A[3,:]
        Nsub = int(N/4.)

        
        mean_vi = np.mean(vi2)
        sd_vi = np.std(vi2)

        maxXvi = mean_vi + 3*sd_vi
        minXvi = mean_vi - 2.5*sd_vi

        Xvi = np.linspace(minXvi, maxXvi, 200)
        # [a, c, scale, BIC]  = estim_param_init_adder(data_or_cdt, author)
        # [k_i, tet_i] = estim_param_init_adder2(data_or_cdt, author)
        Di = 2*vi_2 - vi_1
        Di = Di[Di > 0]
        f_Di = stats.gaussian_kde(Di)
        hatF_Visim = np.zeros((nb_sim, 4, 200))
        for n in range(nb_sim):
            Delta_i = f_Di.resample(N).T.flatten()
            #Delta_i = stats.gamma.rvs(size = N, a = k_i, scale = tet_i)
            Vi_sim = (vi1 + Delta_i)/2.
            hatF_Visim[n, :, :] = funct_rep2(Vi_sim, Nsub, Xvi)
        hatF_Visim = np.sort(hatF_Visim, axis = 0)
        bornesup = hatF_Visim[int((1-alpha/2.)*nb_sim), :, :]
        borneinf = hatF_Visim[int(alpha*nb_sim/2.), :, :]
        F_exp = funct_rep2(vi2, Nsub, Xvi)
        plt.subplot(3, 4, 4*index + 1)
        plt.fill_between(Xvi, borneinf[0,:], bornesup[0,:], color = c[index], alpha = 0.5)
        plt.plot(Xvi, F_exp[0,:], 'k')
        plt.subplot(3, 4, 4*index + 2)
        plt.fill_between(Xvi, borneinf[1,:], bornesup[1,:], color = c[index], alpha = 0.5)
        plt.plot(Xvi, F_exp[1,:], 'k')
        plt.subplot(3, 4, 4*index + 3)
        plt.fill_between(Xvi, borneinf[2,:], bornesup[2,:], color = c[index], alpha = 0.5)
        plt.plot(Xvi, F_exp[2,:], 'k')
        plt.subplot(3, 4, 4*index + 4)
        plt.fill_between(Xvi, borneinf[3,:], bornesup[3,:], color = c[index], alpha = 0.5)
        plt.plot(Xvi, F_exp[3,:], 'k')
        

def figure_simu(DATA, AUTHOR, nb_sim, alpha, index):
    # DATA = [data_Tiruvadi_STK13_Ace, data_Tiruvadi_STK13_AlaTrE, data_si_mops_acet, data_Tiruvadi_JM85_Gly, data_Tiruvadi_STK13_Gly, data_Tiruvadi_STK13_GlyTrE, data_Tiruvadi_STK13_Glu, data_si_mops_arg, 'glycerol', data_Tiruvadi_STK13_GlyCas, data_si_mops_gly11aa, 'glucose', data_si_mops_glu1, 'glucose8a',  data_si_mops_glu2, data_si_mops_glu12aa]
    # AUTHOR = ['tiru', 'tiru', 'si', 'tiru', 'tiru', 'tiru', 'tiru', 'si', 'witz', 'tiru', 'si', 'witz', 'si', 'witz', 'si', 'si']
    figure_simu_DIAM20(DATA, AUTHOR, nb_sim, alpha, index)
    figure_simu_CPM20(DATA, AUTHOR, nb_sim, alpha, index)
    figure_simu_RDAM20(DATA, AUTHOR, nb_sim, alpha, index)
    figure_simu_CAM20(DATA, AUTHOR, nb_sim, alpha, index)
    
    

def fig_simu_iniadder(DATA, AUTHOR, nb_sim, alpha):
    figure_simu_initadder(DATA, AUTHOR, nb_sim, alpha, 0)
    figure_simu_initadder(DATA, AUTHOR, nb_sim, alpha, 1)
    figure_simu_initadder(DATA, AUTHOR, nb_sim, alpha, 2)
    
def Wasserstein_RDAM(data_or_cdt, author, nb_sim):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    N = np.size(v)
    v_tild = np.reshape(v_tild, (N,1))
    V = np.reshape(V, (N,1))
    v = np.reshape(v, (N,1))
    lamb = np.reshape(lamb, (N,1))
    dist_vi = []
    dist_vb = []
    dist_lamb = []
    [a, c, scale, _]  = estim_param_RDAM(data_or_cdt, author)
    for n in range(nb_sim):
        Delta_id = stats.gengamma.rvs(size = (N,1), a = a, c = c, scale = scale)
        Vb_sim = v_tild + Delta_id
        dist_vi.append(wasserstein_distance(np.concatenate((Vb_sim, v_tild), axis = 1), np.concatenate((V, v_tild), axis = 1)))
        dist_vb.append(wasserstein_distance(np.concatenate((Vb_sim, v), axis = 1), np.concatenate((V, v), axis = 1)))
        dist_lamb.append(wasserstein_distance(np.concatenate((Vb_sim, lamb), axis = 1), np.concatenate((V, lamb), axis = 1)))
    dist_vi = np.array(dist_vi)
    dist_vb = np.array(dist_vb)
    dist_lamb = np.array(dist_lamb)
    return np.array([[np.mean(dist_vi), np.std(dist_vi)], 
                     [np.mean(dist_vb), np.std(dist_vb)], 
                     [np.mean(dist_lamb), np.std(dist_lamb)]])

def Wasserstein_DIAM(data_or_cdt, author, nb_sim):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    N = np.size(v)
    v_tild = np.reshape(v_tild, (N,1))
    V = np.reshape(V, (N,1))
    v = np.reshape(v, (N,1))
    lamb = np.reshape(lamb, (N,1))
    dist_vi = []
    dist_vb = []
    dist_lamb = []
    [a, c, scale, _]  = estim_param_DIAM(data_or_cdt, author)
    for n in range(nb_sim):
        Delta_d = stats.gengamma.rvs(size = (N,1), a = a, c = c, scale = scale)
        Vb_sim = (v + Delta_d)/2.
        dist_vi.append(wasserstein_distance(np.concatenate((Vb_sim, v_tild), axis = 1), np.concatenate((V, v_tild), axis = 1)))
        dist_vb.append(wasserstein_distance(np.concatenate((Vb_sim, v), axis = 1), np.concatenate((V, v), axis = 1)))
        dist_lamb.append(wasserstein_distance(np.concatenate((Vb_sim, lamb), axis = 1), np.concatenate((V, lamb), axis = 1)))
    dist_vi = np.array(dist_vi)
    dist_vb = np.array(dist_vb)
    dist_lamb = np.array(dist_lamb)
    return np.array([[np.mean(dist_vi), np.std(dist_vi)], 
                     [np.mean(dist_vb), np.std(dist_vb)], 
                     [np.mean(dist_lamb), np.std(dist_lamb)]])

def Wasserstein_CPM(data_or_cdt, author, nb_sim):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    N = np.size(v)
    v_tild = np.reshape(v_tild, (N,1))
    V = np.reshape(V, (N,1))
    v = np.reshape(v, (N,1))
    lamb = np.reshape(lamb, (N,1))
    dist_vi = []
    dist_vb = []
    dist_lamb = []
    [a_d, s_d, a_r, s_r]  = estim_param_AND(data_or_cdt, author)[0]
    for n in range(nb_sim):
        R = stats.gamma.rvs(size = (N,1), a = a_r, scale = s_r)
        Delta_d = stats.gamma.rvs(size = (N,1), a = a_d, scale = s_d)
        Rep = v_tild*np.exp(lamb*R)
        Div = (v + Delta_d)/2.
        cond_rep = np.where(Rep > Div, 1, 0)
        Vb_sim = Rep*cond_rep + Div*(1 - cond_rep)
        dist_vi.append(wasserstein_distance(np.concatenate((Vb_sim, v_tild), axis = 1), np.concatenate((V, v_tild), axis = 1)))
        dist_vb.append(wasserstein_distance(np.concatenate((Vb_sim, v), axis = 1), np.concatenate((V, v), axis = 1)))
        dist_lamb.append(wasserstein_distance(np.concatenate((Vb_sim, lamb), axis = 1), np.concatenate((V, lamb), axis = 1)))
    dist_vi = np.array(dist_vi)
    dist_vb = np.array(dist_vb)
    dist_lamb = np.array(dist_lamb)
    # plt.plot(*stats.probplot(dist_vi, dist="norm")[0], 'bo')
    return np.array([[np.mean(dist_vi), np.std(dist_vi)], 
                      [np.mean(dist_vb), np.std(dist_vb)], 
                      [np.mean(dist_lamb), np.std(dist_lamb)]])

def Wasserstein_CAM(data_or_cdt, author, nb_sim):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    N = np.size(v)
    v_tild = np.reshape(v_tild, (N,1))
    V = np.reshape(V, (N,1))
    v = np.reshape(v, (N,1))
    lamb = np.reshape(lamb, (N,1))
    dist_vi = []
    dist_vb = []
    dist_lamb = []
    [a_d, s_d, a_r, s_r]  = estim_param_revolution(data_or_cdt, author)[0]
    for n in range(nb_sim):
        Delta_id = stats.gamma.rvs(size = (N,1), a = a_r, scale = s_r)
        Delta_d = stats.gamma.rvs(size = (N,1), a = a_d, scale = s_d)
        Rep = v_tild + Delta_id
        Div = (v + Delta_d)/2.
        cond_rep = np.where(Rep > Div, 1, 0)
        Vb_sim = Rep*cond_rep + Div*(1 - cond_rep)
        dist_vi.append(wasserstein_distance(np.concatenate((Vb_sim, v_tild), axis = 1), np.concatenate((V, v_tild), axis = 1)))
        dist_vb.append(wasserstein_distance(np.concatenate((Vb_sim, v), axis = 1), np.concatenate((V, v), axis = 1)))
        dist_lamb.append(wasserstein_distance(np.concatenate((Vb_sim, lamb), axis = 1), np.concatenate((V, lamb), axis = 1)))
    dist_vi = np.array(dist_vi)
    dist_vb = np.array(dist_vb)
    dist_lamb = np.array(dist_lamb)
    # plt.plot(*stats.probplot(dist_vi, dist="norm")[0], 'bo')
    return np.array([[np.mean(dist_vi), np.std(dist_vi)], 
                      [np.mean(dist_vb), np.std(dist_vb)], 
                      [np.mean(dist_lamb), np.std(dist_lamb)]])



def Wasserstein(nb_sim):
    DATA = [data_Tiruvadi_STK13_Ace, data_Tiruvadi_STK13_AlaTrE, data_si_mops_acet, data_Tiruvadi_JM85_Gly, data_Tiruvadi_STK13_Gly, data_Tiruvadi_STK13_GlyTrE, data_Tiruvadi_STK13_Glu, data_si_mops_arg, 'glycerol', data_Tiruvadi_STK13_GlyCas, data_si_mops_gly11aa, 'glucose', data_si_mops_glu1, 'glucose8a',  data_si_mops_glu2, data_si_mops_glu12aa]
    AUTHOR = ['tiru', 'tiru', 'si', 'tiru', 'tiru', 'tiru', 'tiru', 'si', 'witz', 'tiru', 'si', 'witz', 'si', 'witz', 'si', 'si']
    for i in range(len(DATA)): 
        np.random.seed(123)
        # wass_DIAM = Wasserstein_DIAM(DATA[i], AUTHOR[i], nb_sim)
        # wass_RDAM = Wasserstein_RDAM(DATA[i], AUTHOR[i], nb_sim)
        # wass_CPM = Wasserstein_CPM(DATA[i], AUTHOR[i], nb_sim)
        # Result_i = np.concatenate((wass_DIAM, wass_RDAM, wass_CPM), axis = 1)*100
        Result_i = Wasserstein_CAM(DATA[i], AUTHOR[i], nb_sim)*100
        print(Result_i)
        if i == 0:
            ResultArray = Result_i
        else :
            ResultArray = np.concatenate((ResultArray, Result_i), axis = 0)
    return ResultArray


Wass = np.array([[6.475701, 0.436681, 2.730112, 0.207377, 3.732797, 0.252243],
       [5.429091, 0.493977, 3.297957, 0.211471, 3.117716, 0.217656],
       [4.587729, 0.580838, 2.009547, 0.282355, 1.905706, 0.326868],
       [8.646025, 0.947648, 3.335293, 0.268519, 4.266575, 0.418169],
       [9.124081, 0.942926, 4.386816, 0.257752, 5.296783, 0.444221],
       [7.876631, 1.101487, 2.869929, 0.29654 , 3.226901, 0.408942],
       [8.045785, 0.300276, 3.893995, 0.225675, 5.129964, 0.242693],
       [6.643737, 0.350292, 5.656212, 0.265038, 4.109881, 0.251387],
       [6.142772, 0.358782, 3.48902 , 0.25057 , 3.887437, 0.268323],
       [6.323218, 0.449166, 3.135374, 0.212667, 4.341336, 0.347592],
       [5.757074, 0.471297, 5.082883, 0.346226, 4.290887, 0.289353],
       [4.590218, 0.535692, 2.237207, 0.245207, 2.88672 , 0.320501],
       [5.186089, 0.476957, 4.150835, 0.361831, 5.008039, 0.529654],
       [5.842082, 0.513026, 5.144058, 0.342954, 5.385285, 0.495228],
       [4.51603 , 0.522676, 3.831675, 0.386916, 4.285252, 0.567974],
       [4.899179, 0.472611, 2.780488, 0.311611, 3.460497, 0.395042],
       [5.557518, 0.504672, 4.351623, 0.357665, 4.150139, 0.396192],
       [4.157542, 0.510478, 2.741821, 0.304947, 3.047016, 0.395649],
       [5.911727, 0.863156, 3.808185, 0.40292 , 4.296559, 0.546861],
       [7.166076, 0.745611, 5.145344, 0.380291, 6.017282, 0.587009],
       [5.883482, 0.857335, 4.214934, 0.40847 , 4.1444  , 0.55139 ],
       [6.905322, 0.23922 , 1.561391, 0.114896, 1.718298, 0.134962],
       [6.860942, 0.244381, 2.244302, 0.127396, 1.803698, 0.128172],
       [6.835642, 0.243106, 2.449846, 0.100974, 1.676423, 0.129611],
       [6.948516, 0.385277, 3.084445, 0.276011, 2.849254, 0.214745],
       [6.850597, 0.398208, 4.551916, 0.297368, 3.187388, 0.185905],
       [6.576771, 0.426125, 3.339412, 0.264546, 2.591242, 0.213446],
       [14.438218,  1.133887,  4.915052,  0.612186,  6.743579,  0.832106],
       [14.770715,  1.059867, 14.92473 ,  1.24748 ,  8.903578,  0.730781],
       [15.294829,  0.975505, 12.17054 ,  0.732012,  7.549174,  0.845214],
       [4.904591,  0.361567,  4.734455,  0.451907,  4.008717,  0.370866],
       [5.897251,  0.298934, 13.358719,  0.578469,  5.152409,  0.31192 ],
       [5.700927,  0.33976 ,  6.099289,  0.413374,  5.030376,  0.357733],
       [5.78499 ,  0.470947,  5.308022,  0.427889,  3.882437,  0.333327],
       [6.518938,  0.411566, 10.078197,  0.451309,  5.064042,  0.338384],
       [6.41045 ,  0.426843,  6.069406,  0.384697,  4.143252,  0.318091],
       [4.062483,  0.312862,  3.793153,  0.397399,  3.377951,  0.281603],
       [4.662916,  0.287051, 10.258921,  0.474075,  4.068652,  0.257739],
       [4.660845,  0.334917,  4.814379,  0.400061,  4.472589,  0.294877],
       [9.140783, 0.648602, 4.613051, 0.422813, 4.67134 , 0.442492],
       [9.633358, 0.585289, 9.070166, 0.49098 , 6.56316 , 0.458945],
       [9.722193, 0.589864, 6.71068 , 0.381131, 5.24505 , 0.412703],
       [4.52759 , 0.344045, 1.429159, 0.113017, 1.826733, 0.201799],
       [4.818953, 0.345097, 4.877253, 0.231122, 2.423892, 0.19455 ],
       [4.555643, 0.338955, 2.650235, 0.124653, 1.948762, 0.189224],
       [ 3.528027,  0.500448,  3.168032,  0.39256 ,  1.965957,  0.294616],
       [ 4.587293,  0.388783, 10.97892 ,  0.409945,  3.589724,  0.235186],
       [ 3.667385,  0.477268,  3.698379,  0.344107,  2.180157,  0.286257]])

Wass2 = np.array([[ 4.94158588,  0.32792793,  2.69905073,  0.21010758,  3.74627126, 0.28150585,  2.48569375,  0.17007492],
       [ 3.457736  ,  0.28205851,  3.27279477,  0.19181256,  3.09722837, 0.22855798,  2.96131654,  0.17737696],
       [ 2.13174195,  0.35302736,  1.96849992,  0.27219639,  1.91285583, 0.34638788,  1.62081685,  0.24859154],
       [ 5.16582272,  0.48201352,  3.36821788,  0.22796866,  4.2456812 , 0.41907689,  3.32644333,  0.26997981],
       [ 5.88163195,  0.5520738 ,  4.39588347,  0.24870174,  5.30671592, 0.42029027,  4.33035857,  0.23987529],
       [ 3.63464185,  0.5488759 ,  2.87448319,  0.25228461,  3.21603891, 0.42850594,  2.88306611,  0.27866774],
       [ 6.30141156,  0.20470427,  3.90303438,  0.22332709,  5.12603986, 0.23014176,  2.13644254,  0.19601858],
       [ 4.61937106,  0.23782791,  5.58801928,  0.29245968,  4.13335012, 0.24626172,  2.86589024,  0.17646835],
       [ 4.02158498,  0.23694089,  3.50011085,  0.27231103,  3.8958671 , 0.2598459 ,  1.84830658,  0.17460171],
       [ 5.19284547,  0.3650296 ,  3.1449635 ,  0.1905544 ,  4.36878334, 0.31312504,  3.13514792,  0.19188367],
       [ 4.60250211,  0.35741116,  5.15524826,  0.3550143 ,  4.24592772, 0.30498649,  4.05934805,  0.20945218],
       [ 3.0472667 ,  0.38569108,  2.26621897,  0.26068894,  2.84932081, 0.33067385,  2.22897717,  0.2110308 ],
       [ 4.65751951,  0.38321222,  4.20716818,  0.37019293,  4.9867461 , 0.57311659,  3.90003786,  0.34618056],
       [ 5.00138233,  0.41842288,  5.18910949,  0.34711444,  5.41928157, 0.51652271,  4.41372951,  0.31865254],
       [ 3.79181637,  0.43298709,  3.90211682,  0.38962425,  4.28840033, 0.6265039 ,  3.36127934,  0.38176892],
       [ 4.09032106,  0.3754126 ,  2.81889289,  0.30463328,  3.45034656, 0.37534752,  2.11411235,  0.14347612],
       [ 4.62589529,  0.40801895,  4.38266223,  0.33701325,  4.1556081 , 0.36753961,  2.91613421,  0.13161859],
       [ 3.23347984,  0.39742627,  2.79624616,  0.28840602,  3.03194264, 0.37469759,  2.09029975,  0.15882885],
       [ 4.07609732,  0.45871886,  3.74943562,  0.4505627 ,  4.33592497, 0.51537674,  3.43939726,  0.33764356],
       [ 5.72670151,  0.51500902,  5.13730941,  0.41155301,  6.04137771, 0.53818653,  4.70630757,  0.31125083],
       [ 4.07908414,  0.44963766,  4.18670488,  0.42785596,  4.19537419, 0.51592788,  3.70911476,  0.34656364],
       [ 2.42868564,  0.11329115,  1.59847017,  0.15330743,  1.71552734, 0.13097919,  1.1065975 ,  0.05615321],
       [ 1.49113088,  0.10931347,  2.26043578,  0.13938409,  1.7961672 , 0.12030867,  1.32394812,  0.05626381],
       [ 1.62229556,  0.10373211,  2.47990844,  0.12567698,  1.67474505, 0.12640941,  1.91922264,  0.07695371],
       [ 4.13683162,  0.23459793,  3.10192859,  0.24798289,  2.85635264, 0.22387916,  2.32405564,  0.14113603],
       [ 3.36464625,  0.19036171,  4.61025181,  0.26733369,  3.22760007, 0.1961738 ,  2.85211663,  0.1320971 ],
       [ 2.88064258,  0.20265699,  3.35364054,  0.24527566,  2.613237  , 0.20731323,  2.74651228,  0.16183533],
       [ 5.2814492 ,  0.74869456,  4.88464506,  0.58925997,  6.59688912, 0.86974611,  4.90682871,  0.69374408],
       [ 7.66031042,  0.72570385, 14.8098338 ,  1.27642041,  8.80083374, 0.80826339,  7.55204214,  0.78544436],
       [ 8.7158069 ,  0.81734641, 12.12930262,  0.75987266,  7.43362032, 0.86323065,  8.64260705,  0.80564612],
       [ 4.27741202,  0.35433601,  4.72357171,  0.45223012,  4.02861021, 0.34536   ,  2.89743193,  0.23765288],
       [ 5.30195577,  0.30171738, 13.40491972,  0.59552527,  5.14315416, 0.29432675,  4.86729779,  0.26798987],
       [ 5.0857174 ,  0.33953519,  6.07818326,  0.3929925 ,  5.06375094, 0.32342844,  4.05671501,  0.24181825],
       [ 4.013938  ,  0.33563948,  5.34858125,  0.43044338,  3.94438   , 0.36510915,  3.72210974,  0.30162734],
       [ 5.29745482,  0.33838755, 10.09572273,  0.46215899,  5.10832759, 0.35016874,  4.88208462,  0.30275506],
       [ 4.94517731,  0.31974226,  6.09448341,  0.41980028,  4.2125499 , 0.35747536,  4.80445464,  0.30422485],
       [ 3.588239  ,  0.31720441,  3.8011965 ,  0.35112414,  3.35594886, 0.32346029,  2.51897466,  0.19506818],
       [ 4.14055043,  0.28410167, 10.31150298,  0.37824616,  4.04158287, 0.28148947,  3.22191366,  0.16349959],
       [ 4.27901089,  0.32294964,  4.82356037,  0.35128582,  4.47376533, 0.32412891,  3.24102265,  0.22400249],
       [ 4.60172272,  0.39621434,  4.65097232,  0.44003081,  4.68252438, 0.46447631,  3.39300154,  0.31916727],
       [ 6.09679467,  0.39977817,  9.05112039,  0.48399042,  6.51813633, 0.46404128,  4.64960504,  0.30722631],
       [ 5.87457403,  0.36279293,  6.69509377,  0.36844743,  5.23911493, 0.45604783,  5.32925385,  0.29218493],
       [ 1.94956723,  0.21686369,  1.4192955 ,  0.11179203,  1.8137207 , 0.19325899,  1.45275962,  0.11448305],
       [ 2.44141129,  0.21963042,  4.91862188,  0.23670062,  2.40545008, 0.19149545,  2.17991271,  0.13899841],
       [ 2.01072265,  0.20848094,  2.61723187,  0.14394089,  1.93600488, 0.18643316,  1.89582041,  0.11570961],
       [ 2.7390708 ,  0.30698111,  3.16752938,  0.39141011,  1.98349294, 0.27626793,  2.16560997,  0.24839924],
       [ 4.06774957,  0.22224082, 10.94997839,  0.43288562,  3.57475671, 0.23439816,  4.0106061 ,  0.22169941],
       [ 2.90410034,  0.2876308 ,  3.70786383,  0.34787074,  2.19694605, 0.2717706 ,  2.44375543,  0.23253614]])


def plot_wass():
    Wass3 = Wass2/100.
    t_interval = stats.t.ppf(0.975, df = 99)*np.sqrt(101/100)
    Vb_vi_mean = Wass3[::3,::2]
    Vb_vi_mean = (Vb_vi_mean - np.mean(Vb_vi_mean, axis = 1)[:,None])/np.std(Vb_vi_mean, axis = 1)[:,None]
    Vb_vi_std = Wass3[::3,1::2]*t_interval
    Vb_vb_mean = Wass3[1::3,::2]
    Vb_vb_mean = (Vb_vb_mean - np.mean(Vb_vb_mean, axis = 1)[:,None])/np.std(Vb_vb_mean, axis = 1)[:,None]
    Vb_vb_std = Wass3[1::3,1::2]*t_interval
    Vb_lamb_mean = Wass3[2::3,::2]
    Vb_lamb_std = Wass3[2::3,1::2]*t_interval
    #Lamb = [0.08, 0.21, 0.23, 0.26, 0.27, 0.30, 0.38, 0.41, 0.44, 0.51, 0.57, 0.7, 0.82, 0.92, 1, 1.46]
    XDIAM = [1.1*i - 0.2 for i in range(1,4)] + [1.1*i - 0.2 +1 for i in range(4,7)] + [1.1*i - 0.2 + 2 for i in range(7,10)] + [1.1*i - 0.2 + 3 for i in range(10,17)]
    XCPM = [1.1*i for i in range(1,4)] + [1.1*i+1 for i in range(4,7)] + [1.1*i + 2 for i in range(7,10)] + [1.1*i + 3 for i in range(10,17)]
    XRDAM = [1.1*i + 0.2 for i in range(1,4)] + [1.1*i + 0.2 +1 for i in range(4,7)] + [1.1*i + 0.2 + 2 for i in range(7,10)] + [1.1*i + 0.2 + 3 for i in range(10,17)]
    XCAM = [1.1*i + 0.4 for i in range(1,4)] + [1.1*i + 0.4 +1 for i in range(4,7)] + [1.1*i + 0.4 + 2 for i in range(7,10)] + [1.1*i + 0.4 + 3 for i in range(10,17)]
    strain = np.array(['stk', 'stk', 'mg', 'jm', 'stk', 'stk', 'stk', 'ncm', 'bw', 'stk', 'mg', 'bw', 'mg', 'bw', 'ncm', 'ncm'])
    index_strain = [j for j in range(len(strain)) if strain[j] == 'mg'] + [j for j in range(len(strain)) if strain[j] == 'ncm'] + [j for j in range(len(strain)) if strain[j] == 'bw'] + [j for j in range(len(strain)) if strain[j] == 'stk'] + [j for j in range(len(strain)) if strain[j] == 'jm'] 
    # plt.subplot(311)
    # plt.yscale('log')
    # plt.errorbar(XDIAM, Vb_vi_mean[index_strain,0], yerr = Vb_vi_std[index_strain,0], ecolor = 'k', color = 'palevioletred', capsize=2,fmt='s', alpha = 1, markersize=4,  label = 'DIAM')
    # plt.errorbar(XRDAM, Vb_vi_mean[index_strain,1], yerr = Vb_vi_std[index_strain,1], ecolor = 'k', color = 'darkcyan', capsize=2,fmt='s', alpha = 1, markersize=4,  label = 'RDAM')
    # plt.errorbar(XCPM, Vb_vi_mean[index_strain,2], yerr = Vb_vi_std[index_strain,2], ecolor = 'k', color = 'slateblue', capsize=2,fmt='s', alpha = 1, markersize=4, label = 'CPM')
    # plt.errorbar(XCAM, Vb_vi_mean[index_strain,3], yerr = Vb_vi_std[index_strain,3], ecolor = 'k', color = 'darkred', capsize=2,fmt='s', alpha = 1, markersize=4, label = 'CAM')
    # plt.ylim(0.009, 0.2)
    # plt.subplot(312)
    # plt.yscale('log')
    # plt.errorbar(XDIAM, Vb_vb_mean[index_strain,0], yerr = Vb_vb_std[index_strain,0], ecolor = 'k', color = 'palevioletred', capsize=2,fmt='s', alpha = 1, markersize=4,  label = 'DIAM')
    # plt.errorbar(XRDAM, Vb_vb_mean[index_strain,1], yerr = Vb_vb_std[index_strain,1], ecolor = 'k', color = 'darkcyan', capsize=2,fmt='s', alpha = 1, markersize=4,  label = 'RDAM')
    # plt.errorbar(XCPM, Vb_vb_mean[index_strain,2], yerr = Vb_vb_std[index_strain,2], ecolor = 'k', color = 'slateblue', capsize=2,fmt='s', alpha = 1, markersize=4, label = 'CPM')
    # plt.errorbar(XCAM, Vb_vb_mean[index_strain,3], yerr = Vb_vb_std[index_strain,3], ecolor = 'k', color = 'darkred', capsize=2,fmt='s', alpha = 1, markersize=4, label = 'CAM')
    # plt.yscale('log')
    # plt.errorbar(XDIAM, Vb_lamb_mean[index_strain,0], yerr = Vb_lamb_std[index_strain,0], ecolor = 'k', color = 'palevioletred', capsize=2,fmt='s', alpha = 1, markersize=4,  label = 'DIAM')
    # plt.errorbar(XRDAM, Vb_lamb_mean[index_strain,1], yerr = Vb_lamb_std[index_strain,1], ecolor = 'k', color = 'darkcyan', capsize=2,fmt='s', alpha = 1, markersize=4,  label = 'RDAM')
    # plt.errorbar(XCPM, Vb_lamb_mean[index_strain,2], yerr = Vb_lamb_std[index_strain,2], ecolor = 'k', color = 'slateblue', capsize=2,fmt='s', alpha = 1, markersize=4, label = 'CPM')
    # plt.errorbar(XCAM, Vb_lamb_mean[index_strain,3], yerr = Vb_lamb_std[index_strain,3], ecolor = 'k', color = 'darkred', capsize=2,fmt='s', alpha = 1, markersize=4, label = 'CAM')
    
    plt.scatter(Vb_vi_mean[:,0], Vb_vb_mean[:,0],  color = 'palevioletred', alpha = 1)
    plt.scatter(Vb_vi_mean[:,1], Vb_vb_mean[:,1],  color = 'darkcyan', alpha = 1)
    plt.scatter(Vb_vi_mean[:,2], Vb_vb_mean[:,2],  color = 'slateblue', alpha = 1)
    plt.scatter(Vb_vi_mean[:,3], Vb_vb_mean[:,3],  color = 'darkred', alpha = 1)
    

def wasserstein_distance(X, Y):
    n_X = np.size(X[:,0])
    n_Y = np.size(Y[:,0])
    cost_matrix = ot.dist(X, Y, metric='euclidean')
    weights_X = np.ones(n_X) / n_X
    weights_Y = np.ones(n_Y) / n_Y
    return ot.emd2(weights_X, weights_Y, cost_matrix)


def test_estim_param_AND(data_or_cdt, author):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    E_Dd = np.mean(Delta_app)
    E_R= np.mean(R_app)
    results = estim_param_AND(data_or_cdt, author)
    print(results)
    [a_d, s_d, a_r, s_r] = results[0]
    A = np.linspace(0, E_Dd*3, 100)
    R = np.linspace(0, E_R*3, 100)
    f_delta_app = []
    for dd in A:
        f_delta_app.append(np.mean(stats.gamma.pdf(dd, a = a_d, scale = s_d)*stats.gamma.cdf(R_app + np.log((v + dd)/(2*V))/lamb, a = a_r, scale = s_r) + stats.gamma.cdf(dd, a = a_d, scale = s_d)*stats.gamma.pdf(R_app + np.log((v + dd)/(2*V))/lamb, a = a_r, scale = s_r)/(lamb*(v + dd))))
    f_r_app = []
    for rr in R:
        aa = 2*v_tild*np.exp(lamb*rr)
        f_r_app.append(np.mean(aa*lamb*stats.gamma.pdf(aa - v, a = a_d, scale = s_d)*stats.gamma.cdf(rr, a = a_r, scale = s_r) + stats.gamma.cdf(aa - v, a = a_d, scale = s_d)*stats.gamma.pdf(rr, a = a_r, scale = s_r)))
    plt.subplot(211)
    plt.hist(Delta_app, [3*k*E_Dd/50. for k in range(50)], color = 'k', density = True)
    plt.plot(A, f_delta_app, 'c', label = '$\Delta_{app}$')
    plt.plot(A, stats.gamma.pdf(A, a = a_d, scale = s_d), 'r', label = '$\Delta_d$')
    plt.legend()
    plt.ylabel('$Density$')
    plt.xlabel('$Added$ $size$ ($\mu m$)')
    plt.subplot(212)
    plt.hist(R_app, [3*k*E_R/50. for k in range(50)], color = 'k', density = True)
    plt.plot(R, f_r_app, 'c', label = '$R_{app}$')
    plt.plot(R, stats.gamma.pdf(R, a = a_r, scale = s_r), 'r', label = '$R$')
    plt.legend()
    plt.ylabel('$Density$')
    plt.xlabel('$C+D$ $period$ ($h$)')

def test_estim_revolution(data_or_cdt, author):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    E_Dd = np.mean(Delta_app)
    E_Did= np.mean(V - v_tild)
    results = estim_param_revolution(data_or_cdt, author)
    print(results)
    [a_d, s_d, a_r, s_r] = results[0]
    A = np.linspace(0, E_Dd*3, 100)
    R = np.linspace(0, E_Did*3, 100)
    f_delta_app = []
    for dd in A:
        f_delta_app.append(np.mean(stats.gamma.pdf(dd, a = a_d, scale = s_d)*stats.gamma.cdf((dd + v)/2. - v_tild, a = a_r, scale = s_r) + stats.gamma.cdf(dd, a = a_d, scale = s_d)*stats.gamma.pdf((dd + v)/2. - v_tild, a = a_r, scale = s_r)/2.))
    f_r_app = []
    for rr in R:
        aa = 2*(rr + v_tild)
        f_r_app.append(np.mean(2*stats.gamma.pdf(aa - v, a = a_d, scale = s_d)*stats.gamma.cdf(rr, a = a_r, scale = s_r) + stats.gamma.cdf(aa - v, a = a_d, scale = s_d)*stats.gamma.pdf(rr, a = a_r, scale = s_r)))
    plt.subplot(211)
    plt.hist(Delta_app, [3*k*E_Dd/50. for k in range(50)], color = 'k', density = True)
    plt.plot(A, f_delta_app, 'c', label = '$\Delta_{app}$')
    plt.plot(A, stats.gamma.pdf(A, a = a_d, scale = s_d), 'r', label = '$\Delta_d$')
    plt.legend()
    plt.ylabel('$Density$')
    plt.xlabel('$Added$ $size$ ($\mu m$)')
    plt.subplot(212)
    plt.hist(V - v_tild, [3*k*E_Did/50. for k in range(50)], color = 'k', density = True)
    plt.plot(R, f_r_app, 'c', label = '$\Delta_{id app}$')
    plt.plot(R, stats.gamma.pdf(R, a = a_r, scale = s_r), 'r', label = '$\Delta_{id}$')
    plt.legend()
    plt.ylabel('$Density$')
    plt.xlabel('$C+D$ $period$ ($h$)')

def estim_param_DIAM(data_or_cdt, author):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    a, c ,_,scale = stats.gengamma.fit(Delta_app, floc = 0)
    BIC = - np.sum(np.log(2*stats.gengamma.pdf(Delta_app, a = a, c = c, scale = scale))) + 3*np.log(np.size(v))/2.
    return [a, c, scale, BIC]

def estim_param_RDAM(data_or_cdt, author):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    a, c ,_,scale = stats.gengamma.fit(V - v_tild, floc = 0)
    BIC = - np.sum(np.log(stats.gengamma.pdf(V - v_tild, a = a, c = c, scale = scale))) + 3*np.log(np.size(v))/2.
    return [a, c, scale, BIC]

def estim_param_init_adder(data_or_cdt, author):
    [_, Delta_i] = extract_donnees3(data_or_cdt, author)
    moy = np.mean(Delta_i)
    var = np.var(Delta_i)
    #a, c ,_,scale = stats.gengamma.fit(Delta_i, floc = 0)
    # BIC = - np.sum(np.log(stats.gengamma.pdf(Delta_i, a = a, c = c, scale = scale))) + 3*np.log(np.size(Delta_i))/2.
    # return [a, c, scale, BIC]
    return [moy, var]


def estim_param_init_adder2(data_or_cdt, author):
    [vb_1, vi_1, vi_2] = extract_donnees4(data_or_cdt, author)
    Delta_i = 2*vi_2 - vi_1
    E_Di = np.mean(Delta_i)
    V_Di = np.var(Delta_i)
    k_i = E_Di**2 / V_Di
    tet_i = V_Di / E_Di
    return [k_i, tet_i]


def calcul_taux_el(data_or_cdt, author):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    return np.mean(np.log(V/v_tild)/R_app)

def calcul_moy_vb(data_or_cdt, author):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    return np.mean(v)

def calcul_moy_vi(data_or_cdt, author):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    return np.moy(v_tild)

def calcul_k_vi(data_or_cdt, author):
    [v, v_tild, V, Delta_app, R_app, lamb] = extract_donnees(data_or_cdt, author)
    return np.mean(v_tild)**2/np.var(v_tild)

def calcul_moy_tau_C(data_or_cdt, author):
    [v, V, v_tild, C, lamb] = extract_donnees2(data_or_cdt, author)
    return np.mean(C)



def adder_init(data):
    cond = np.where(data[1:, data[0] == 'cell ID'] == data[0:-1, data[0] == 'daughter ID'], 1, 0).T[0]
    taille_init_1 = data[np.concatenate((cond, np.array([0])), axis = 0) == 1, data[0] == 'initiation size per ori (micron)']
    taille_init_2 = data[np.concatenate((np.array([0]), cond), axis = 0) == 1, data[0] == 'initiation size per ori (micron)']
    elongation_rate_1 = data[np.concatenate((cond, np.array([0])), axis = 0) == 1, data[0] == 'elongation rate (1/hour)']
    elongation_rate_2 = data[np.concatenate((np.array([0]), cond), axis = 0) == 1, data[0] == 'elongation rate (1/hour)']
    C1 = data[np.concatenate((cond, np.array([0])), axis = 0) == 1, data[0] == 'C period (minute)']
    C2 = data[np.concatenate((np.array([0]), cond), axis = 0) == 1, data[0] == 'C period (minute)']
    Delta_1 = data[np.concatenate((cond, np.array([0])), axis = 0) == 1, data[0] == 'added size (micron)']
    #V_b_2 = data[np.concatenate((np.array([0]), cond), axis = 0) == 1, data[0] == 'birth size (micron)']
    V_b_2 = 'wrong one'
    return [taille_init_1, taille_init_2, Delta_1, elongation_rate_1, elongation_rate_2, V_b_2, C1, C2]


# Function taking data of Si et al. and returning the data of cells with the information of their mother
def select_data_si_seq_vi(data):
    # indentification of the index of cells having as mother the previous cell
    index = np.where(data[1:, data[0] == 'cell ID'] == data[0:-1, data[0] == 'daughter ID'])[0]
    
    # construction of new array with the information of the mothers
    new_N = np.size(index)
    vi_1 = data[index, data[0] == 'initiation size per ori (micron)']
    vb_1 = data[index, data[0] == 'newborn size (micron)']
    C1 = data[index, data[0] == 'C period (minute)']
    CD1 = data[index, data[0] == 'tau_cyc (minute)']
    added_data = np.reshape(np.concatenate((vi_1, vb_1, C1, CD1)), (4, new_N)).T
    
    # concatenation of the information of the mother with the information of the daughter
    variable_names = np.concatenate((data[0], np.array(['previous initiation size per ori (micron)', 'previous newborn size (micron)', 'previous C period (minute)', 'previous tau_cyc (minute)'])))
    variable_names = np.reshape(variable_names, (1, np.size(variable_names)))
    new_data = np.concatenate((variable_names, np.concatenate((data[index+1,:], added_data), axis = 1)), axis = 0)
    return new_data



def bootstrap_AND(data, author):
    param = []
    N = 1000 #number of iterations
    S = np.shape(data)[0]
    cpt = 0
    for sampling in range(N):
        index = np.concatenate((np.array([0]), np.random.choice(np.arange(1, S), size = S-1)))
        data_sample = data[index, :]
        results = estim_param_AND(data_sample, author)
        if results[1]==True:
            param_div = results[0]
            param.append([param_div[0]*param_div[1], 1/np.sqrt(param_div[0]), param_div[2]*param_div[3], 1/np.sqrt(param_div[2])])
        else :
            print('pas réussie :(', results[2])
            cpt += 1
    
    N = N - cpt
    param = np.array(param)
    
    K_d = np.sort(param[:,0])
    Tet_d = np.sort(param[:,1])
    K_r = np.sort(param[:,2])
    Tet_r = np.sort(param[:,3])

    moy_kd = np.mean(K_d)
    std_kd = np.std(K_d)
    
    moy_tetd = np.mean(Tet_d)
    std_tetd = np.std(Tet_d)
    
    moy_kr = np.mean(K_r)
    std_kr = np.std(K_r)
    
    moy_tetr = np.mean(Tet_r)
    std_tetr = np.std(Tet_r)

    plt.subplot(221)
    plt.plot((K_d - moy_kd)/std_kd, np.sort(stats.norm.rvs(size = N)), 'ko', markerfacecolor = 'none', linestyle = 'none')
    plt.subplot(222)
    plt.plot((Tet_d - moy_tetd)/std_tetd, np.sort(stats.norm.rvs(size = N)), 'ko', markerfacecolor = 'none', linestyle = 'none')
    plt.subplot(223)
    plt.plot((K_r - moy_kr)/std_kr, np.sort(stats.norm.rvs(size = N)), 'ko', markerfacecolor = 'none', linestyle = 'none')
    plt.subplot(224)
    plt.plot((Tet_r - moy_tetr)/std_tetr, np.sort(stats.norm.rvs(size = N)), 'ko', markerfacecolor = 'none', linestyle = 'none')
    
    return [stats.t.interval(0.95, N, loc=moy_kd, scale=std_kd), stats.t.interval(0.95, N, loc=moy_tetd, scale=std_tetd), stats.t.interval(0.95, N, loc=moy_kr, scale=std_kr), stats.t.interval(0.95, N, loc=moy_tetr, scale=std_tetr)]
    

def bootstrap_AND_witzdata(condi):
    data = data_witz[data_witz[condition] == condi]
    data = data.dropna(subset = (Lb, dL, Li, Tid))
    data = data.to_numpy()
    data = data[:,:-1].astype(float)
    param = []
    N = 1000 #number of iterations
    S = np.shape(data)[0]
    cpt = 0
    for sampling in range(N):
        index = np.random.choice(np.arange(0, S), size = S)
        data_sample = data[index, :]
        results = estim_param_ANDgate_witz(data_sample)
        param_div = results[0]
        if results[1] == True:
            param.append([param_div[0]*param_div[1], 1/np.sqrt(param_div[0]), param_div[2]*param_div[3], 1/np.sqrt(param_div[2])])
        else :
            print('pas réussit :(', results[2])
            cpt += 1
    
    N = N - cpt
    param = np.array(param)
    
    K_d = np.sort(param[:,0])
    Tet_d = np.sort(param[:,1])
    K_r = np.sort(param[:,2])
    Tet_r = np.sort(param[:,3])

    moy_kd = np.mean(K_d)
    std_kd = np.std(K_d)
    
    moy_tetd = np.mean(Tet_d)
    std_tetd = np.std(Tet_d)
    
    moy_kr = np.mean(K_r)
    std_kr = np.std(K_r)
    
    moy_tetr = np.mean(Tet_r)
    std_tetr = np.std(Tet_r)

    plt.subplot(221)
    plt.plot((K_d - moy_kd)/std_kd, np.sort(stats.norm.rvs(size = N)), 'ko', markerfacecolor = 'none', linestyle = 'none')
    plt.subplot(222)
    plt.plot((Tet_d - moy_tetd)/std_tetd, np.sort(stats.norm.rvs(size = N)), 'ko', markerfacecolor = 'none', linestyle = 'none')
    plt.subplot(223)
    plt.plot((K_r - moy_kr)/std_kr, np.sort(stats.norm.rvs(size = N)), 'ko', markerfacecolor = 'none', linestyle = 'none')
    plt.subplot(224)
    plt.plot((Tet_r - moy_tetr)/std_tetr, np.sort(stats.norm.rvs(size = N)), 'ko', markerfacecolor = 'none', linestyle = 'none')
    
    return [stats.t.interval(0.95, N, loc=moy_kd, scale=std_kd), stats.t.interval(0.95, N, loc=moy_tetd, scale=std_tetd), stats.t.interval(0.95, N, loc=moy_kr, scale=std_kr), stats.t.interval(0.95, N, loc=moy_tetr, scale=std_tetr)]
    

    

    




   




 
    
