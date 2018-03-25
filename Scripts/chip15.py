
import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d
from VG_ConstPulse_IDSAT import IDSAT
from I_V_curves import hist_IDS_VGS, IDS_VGS_stress, IDS_VGS, IDSAT_vs_row
from VG_ConstPulse_horizontal_hist import IDSAT_horizontal_hist
from VG_ConstPulse_separation import IDSAT_separation
from Charge_Pumping import Charge_Pumping, Charge_Pumping_compare
from MLC_IDSAT import *
from Marcov_Chain import *
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg as sp_lin

def main():
    
    """
    Istep = 0.5
    Imax = 160.0
    Ith = np.array([95.0, 80.0, 65.0, 50.0, 35.0, 20.0, 5.0])
    Tpulse = np.array([0.001, 0.003, 0.009, 0.027, 0.081, 0.243, 0.729])
    Imin = 0.0
    #generator = death_birth_generator_matrix(1/14.0, 40000.0, 0.7, 350, Istep)
    #generator = death_birth_generator_matrix(1/14.0, 100000.0, 1.5, 350, Istep)
    #generator = death_birth_generator_matrix(1/12.0, 120000.0, 5.6, 350, Istep)
    N_state_total = int(round((Imax - Imin)/Istep) + 1)
    k_inv = 15.489
    b0 = 21207.8
    #d0 = 80000.0
    #r = -1.5
    #generator = death_birth_generator_matrix(1/k_inv, b0, d0, r, N_state_total, Istep)

    #generator = Testing_death_birth_generator_matrix(1/k_inv, b0, 4000.0, 13.0, N_state_total, Istep)
    generator = Testing_death_birth_generator_matrix(1/k_inv, b0, 6000.0, 14.0, N_state_total, Istep)

    Prior = Gaussian_function(np.arange(N_state_total), IDSAT_to_state(144.685), 5.778678/Istep)
    relaxation = Gaussian_function(np.arange(N_state_total), IDSAT_to_state(1.35, Imax=160.0, Istep=Istep), 0.95/Istep)
    for cycle in np.arange(6):
        N_state_transient = int(round((Imax - Ith[cycle])/Istep) + 1)
        T_mat, P_mat = transition_matrix(generator, Tpulse[cycle], N_state_total, N_state_transient) 
        Final_mat = np.linalg.matrix_power(P_mat, 1024)
        Prior = np.matmul(Prior, Final_mat)
        Prior_mean = 0
        Prior_std = 0
        for i in np.arange(N_state_total):
            Prior_mean += (Imax - Istep*i) * Prior[i]
        for i in np.arange(N_state_total):
            Prior_std += (Imax - Istep*i - Prior_mean)**2 * Prior[i]
        Prior_std = np.sqrt(Prior_std)
        print('level ', cycle+1, 'mean and std')
        print('before relax', Prior_mean, Prior_std)
        convolve_relax = np.convolve(Prior, relaxation)
        convolve_relax = convolve_relax[-N_state_total:]
        #relax_mean = Imax - Istep * np.average(np.arange(N_state_total), weights = convolve_relax)
        #print(relax_mean)
        relax_mean = 0
        relax_std = 0
        for i in np.arange(N_state_total):
            relax_mean += (Imax - Istep*i) * convolve_relax[i]
        for i in np.arange(N_state_total):
            relax_std += (Imax - Istep*i - relax_mean)**2 * convolve_relax[i]
        relax_std = np.sqrt(relax_std)
        print('after relax', relax_mean, relax_std)
        #plt.figure(cycle+1)
        plt.plot(np.arange(Imax, Imin-1e-4, -Istep), Prior, 'b')
        plt.plot(np.arange(Imax, Imin-1e-4, -Istep), convolve_relax[-N_state_total:], 'r')
        plt.grid()
    """

    Marcov_Chain_MLE(15, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, range(0, 128), [227, 134, 93, 101, 65, 97], [0.001, 0.003, 0.009, 0.027, 0.081, 0.243], 6, ['../Data/chip15/MLC_programming_Chip15_Col33_1msPULSE_VG1p8_VD2p0_VAsource_VBdrain_01', '../Data/chip15/MLC_programming_Chip15_Col33_3msPULSE_VG1p8_VD2p0_VAsource_VBdrain_02', '../Data/chip15/MLC_programming_Chip15_Col33_9msPULSE_VG1p8_VD2p0_VAsource_VBdrain_03', '../Data/chip15/MLC_programming_Chip15_Col33_27msPULSE_VG1p8_VD2p0_VAsource_VBdrain_04', '../Data/chip15/MLC_programming_Chip15_Col33_81msPULSE_VG1p8_VD2p0_VAsource_VBdrain_05', '../Data/chip15/MLC_programming_Chip15_Col33_243msPULSE_VG1p8_VD2p0_VAsource_VBdrain_06'], '../Plots/chip15/', 'VG1p8_VD2p0', '', 160.0, [95.0, 80.0, 65.0, 50.0, 35.0, 20.0], 0.0, 0.5)
    plt.show()

    """
    #fig = plt.figure(1)
    #ax = fig.gca(projection='3d')
    #x = np.arange(350)
    #y = np.arange(350)
    #x, y = np.meshgrid(x, y)
    #surf = ax.plot_surface(x, y, generator)
    ##plt.show()

    #T_mat, P_mat = transition_matrix(generator, 0.01, 350, 325) 
    ##T_mat, P_mat = transition_matrix(generator, 0.001, 350, 325) 
    Final_mat = np.linalg.matrix_power(P_mat, 227)
    #print(np.sum(Final_mat, axis=1))
    #fig = plt.figure(2)
    #ax = fig.gca(projection='3d')
    x = np.arange(350)
    y = np.arange(350)
    x, y = np.meshgrid(x, y)
    #surf = ax.plot_surface(x, y, Final_mat)

    Prior = Gaussian_function(np.arange(350), IDSAT_to_state(144.685), 5.778678/Istep)
    Final_dis = np.matmul(Prior, Final_mat)

    #plt.plot(np.arange(350), Prior)
    #plt.plot(np.arange(350), Final_dis)
    plt.plot(np.arange(Imax, Imin, -Istep), Final_dis)
    plt.xlim(90, 95)
    plt.xlabel('IDSAT (uA)')
    plt.ylabel('distribution')
    #plt.xlim(320, 350)
    plt.grid()
    #print(Final_dis)
    #mean = np.mean(Final_dis)
    #std_dev = np.std(Final_dis)
    #print(mean, std_dev)
    print(np.sum(Final_dis))


    #s_samp = np.zeros((len(np.arange(150.0, 90.0, -1.0))), dtype = np.int32)
    #for (i, I) in zip(np.arange(len(np.arange(150.0, 90.0, -1.0))), np.arange(150.0, 90.0, -1.0)):
    #    s_samp[i] = IDSAT_to_state(I)
    s_samp = np.arange(50, 325, 5) 
    Current = Imax - Istep * s_samp
    #print(s_samp)

    #t_1st = np.zeros((len(s_samp), N_state_total)) 
    #for (i, s) in zip(np.arange(len(s_samp)), s_samp):
    #    T_new = T_mat
    #    T_new[int(s)] = np.zeros(N_state_total)
    #    T_new[int(s)][int(s)] = 1
    #    b = np.ones(N_state_total)
    #    b[int(s)] = 0
    #    t_1st[i] = sp_lin.solve(np.identity(N_state_total) - T_new, b)
    #    print(t_1st[i])

    #t_1st = np.zeros((len(s_samp), N_state_total - 1)) 
    t_1st = np.zeros((len(s_samp), N_state_total)) 
    #for i in np.arange(N_state_total):
    #    print(i, T_mat[i][i])
    for (i, s) in zip(np.arange(len(s_samp)), s_samp):
        t_1st[i][:s] = sp_lin.solve(np.identity(s) - T_mat[:s, :s], np.ones(s))
        #print(t_1st[i][0])
        #print(t_1st[i][50])

    popt, pcov = curve_fit(log_shift, t_1st[:, 50]*Tpulse, Current, bounds=(0, np.inf))
    #popt, pcov = curve_fit(log_shift, t_1st[:, 0]*Tpulse, Current, bounds=(0, np.inf))
    print(popt)
    #rmse
    plt.figure(2)
    plt.plot(t_1st[:, 50]*Tpulse, log_shift(t_1st[:, 50]*Tpulse, *popt), color='m')
    plt.plot(t_1st[:, 50]*Tpulse, Current, marker='.', color = 'b')
    plt.xlabel('programming time (s)')
    plt.ylabel('average IDSAT (uA)')
    #plt.plot(t_1st[:, 0]*Tpulse, log_shift(t_1st[:, 0]*Tpulse, *popt), color='r')
    #plt.plot(t_1st[:, 0]*Tpulse, Current, marker='.', color = 'b')
    plt.grid()
    #plt.figure(2)
    #plt.semilogx(t_1st[:, 0]*Tpulse, Current)
    #plt.grid()

    #for (i, s) in zip(np.arange(len(s_samp)), s_samp):
    #    T_new = np.delete(np.delete(T_mat, s, axis=1), s, axis=0)
    #    t_1st[i] = sp_lin.solve(np.identity(N_state_total - 1) - T_new, np.ones(N_state_total-1))
    #    print(t_1st[i])


    #w, vl, vr = sp_lin.eig(P_mat, left = True, right = True)
    #for i in np.arange(len(w)):
    #    print(i, w[i])
    #for i in np.arange(324, 350):
    #    print('vl '+str(i)+': ', vl[i])
    #    print('vr '+str(i)+': ', vr[i])
    """

    """
    p_b_tot = np.zeros((N_state_transient))
    p_d_tot = np.zeros((N_state_transient))
    p_self = np.zeros((N_state_transient))
    for s in np.arange(N_state_transient):
        p_b_tot[s] = np.sum(P_mat[s][s+1:])
        p_d_tot[s] = np.sum(P_mat[s][:s])
        p_self[s] = P_mat[s][s]

    print(len(p_b_tot), len(np.arange(Imax, Ith-0.01, -Istep)))
    pb, = plt.plot(np.arange(Imax, Ith, -Istep), p_b_tot, marker = '.', color = 'g')
    pd, = plt.plot(np.arange(Imax, Ith, -Istep), p_d_tot, marker = '.', color = 'r')
    pself, = plt.plot(np.arange(Imax, Ith, -Istep), p_self, marker = '.', color = 'y')
    plt.legend([pb, pd, pself], ['total birth', 'total death', 'self transition'], loc='upper left')
    plt.xlabel('current state: IDSAT (uA)')
    plt.ylabel('probability')
    plt.title('L=16, Nfin=2, ULVT, level 1, VG=1.8, VD=2.0, pulse_width 0.1sec\nContinuous Markov model: self transition and total birth and death rate at each transient current state', fontsize=7)
    #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_level-'+str(cycle)+'_P_self_tot-d-b'+title+'.pdf')
    #plt.show()
    """


    #Marcov_Chain_MLE(15, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, range(0, 128), [227], [0.001], 1, ['../Data/chip15/MLC_programming_Chip15_Col33_1msPULSE_VG1p8_VD2p0_VAsource_VBdrain_01'], '../Plots/chip15/', 'VG1p8_VD2p0', '', 160.0, 95.0, 90.0, 0.2)
    #Marcov_Chain_MLE(15, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, range(0, 128), [227], [0.001], 1, ['../Data/chip15/MLC_programming_Chip15_Col33_1msPULSE_VG1p8_VD2p0_VAsource_VBdrain_01'], '../Plots/chip15/', 'VG1p8_VD2p0', '', 160.0, 95.0, 90.0, 1.0)
    #Marcov_Chain_MLE(15, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, range(0, 128), [227], [0.001], 1, ['../Data/chip15/MLC_programming_Chip15_Col33_1msPULSE_VG1p8_VD2p0_VAsource_VBdrain_01'], '../Plots/chip15/', 'VG1p8_VD2p0', '', 160.0, [95.0], 0.0, 0.5)
    #Marcov_Chain_MLE(15, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, range(0, 128), [227, 134, 93, 101, 65, 97, 277], [0.001, 0.003, 0.009, 0.027, 0.081, 0.243, 0.729], 7, ['../Data/chip15/MLC_programming_Chip15_Col33_1msPULSE_VG1p8_VD2p0_VAsource_VBdrain_01', '../Data/chip15/MLC_programming_Chip15_Col33_3msPULSE_VG1p8_VD2p0_VAsource_VBdrain_02', '../Data/chip15/MLC_programming_Chip15_Col33_9msPULSE_VG1p8_VD2p0_VAsource_VBdrain_03', '../Data/chip15/MLC_programming_Chip15_Col33_27msPULSE_VG1p8_VD2p0_VAsource_VBdrain_04', '../Data/chip15/MLC_programming_Chip15_Col33_81msPULSE_VG1p8_VD2p0_VAsource_VBdrain_05', '../Data/chip15/MLC_programming_Chip15_Col33_243msPULSE_VG1p8_VD2p0_VAsource_VBdrain_06', '../Data/chip15/MLC_programming_Chip15_Col33_729msPULSE_VG1p8_VD2p0_VAsource_VBdrain_07'], '../Plots/chip15/', 'VG1p8_VD2p0', '', 160.0, [95.0, 80.0, 65.0, 50.0, 35.0, 20.0, 5.0], 0.0, 0.5)

    #Marcov_Chain(15, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, range(0, 128), [227], [0.001], 1, ['../Data/chip15/MLC_programming_Chip15_Col33_1msPULSE_VG1p8_VD2p0_VAsource_VBdrain_01'], '../Plots/chip15/', 'VG1p8_VD2p0', '', 160.0, 95.0, 90.0, 0.2)

    #for row_start in np.arange(0, 128):
    #    MLC_IDSAT_algorithm_rv1(15, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, [row_start], [227, 134, 93, 101, 65, 97, 277], [0.001, 0.003, 0.009, 0.027, 0.081, 0.243, 0.729], 7, [], '', ['../Data/chip15/MLC_programming_Chip15_Col33_1msPULSE_VG1p8_VD2p0_VAsource_VBdrain_01', '../Data/chip15/MLC_programming_Chip15_Col33_3msPULSE_VG1p8_VD2p0_VAsource_VBdrain_02', '../Data/chip15/MLC_programming_Chip15_Col33_9msPULSE_VG1p8_VD2p0_VAsource_VBdrain_03', '../Data/chip15/MLC_programming_Chip15_Col33_27msPULSE_VG1p8_VD2p0_VAsource_VBdrain_04', '../Data/chip15/MLC_programming_Chip15_Col33_81msPULSE_VG1p8_VD2p0_VAsource_VBdrain_05', '../Data/chip15/MLC_programming_Chip15_Col33_243msPULSE_VG1p8_VD2p0_VAsource_VBdrain_06', '../Data/chip15/MLC_programming_Chip15_Col33_729msPULSE_VG1p8_VD2p0_VAsource_VBdrain_07'], '../Plots/chip15/', 'VG1p8_VD2p0_true-I0_fresh-0p1ms', '_cycle1-7_row'+str(row_start).zfill(3), Imin=0, Imax=160)
    #MLC_IDSAT_algorithm_rv1(15, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, range(0, 128), [227, 134, 93, 101, 65, 97, 277], [0.001, 0.003, 0.009, 0.027, 0.081, 0.243, 0.729], 7, [], '', ['../Data/chip15/MLC_programming_Chip15_Col33_1msPULSE_VG1p8_VD2p0_VAsource_VBdrain_01', '../Data/chip15/MLC_programming_Chip15_Col33_3msPULSE_VG1p8_VD2p0_VAsource_VBdrain_02', '../Data/chip15/MLC_programming_Chip15_Col33_9msPULSE_VG1p8_VD2p0_VAsource_VBdrain_03', '../Data/chip15/MLC_programming_Chip15_Col33_27msPULSE_VG1p8_VD2p0_VAsource_VBdrain_04', '../Data/chip15/MLC_programming_Chip15_Col33_81msPULSE_VG1p8_VD2p0_VAsource_VBdrain_05', '../Data/chip15/MLC_programming_Chip15_Col33_243msPULSE_VG1p8_VD2p0_VAsource_VBdrain_06', '../Data/chip15/MLC_programming_Chip15_Col33_729msPULSE_VG1p8_VD2p0_VAsource_VBdrain_07'], '../Plots/chip15/', 'VG1p8_VD2p0', '_cycle1-7_all')

    #IDS_VGS(15, 33, 16, 2, 'ULVT', 128, ['../Data/chip15/Fresh_Chip15_Col33_Ids_Vgs_VAsource_VBdrain', '../Data/chip15/MLC_Chip15_Col33_1msPULSE_VG1p8_VD2p0_Ids_Vgs_VAsource_VBdrain_01', '../Data/chip15/MLC_Chip15_Col33_3msPULSE_VG1p8_VD2p0_Ids_Vgs_VAsource_VBdrain_02', '../Data/chip15/MLC_Chip15_Col33_9msPULSE_VG1p8_VD2p0_Ids_Vgs_VAsource_VBdrain_03', '../Data/chip15/MLC_Chip15_Col33_27msPULSE_VG1p8_VD2p0_Ids_Vgs_VAsource_VBdrain_04', '../Data/chip15/MLC_Chip15_Col33_81msPULSE_VG1p8_VD2p0_Ids_Vgs_VAsource_VBdrain_05', '../Data/chip15/MLC_Chip15_Col33_243msPULSE_VG1p8_VD2p0_Ids_Vgs_VAsource_VBdrain_06', '../Data/chip15/MLC_Chip15_Col33_729msPULSE_VG1p8_VD2p0_Ids_Vgs_VAsource_VBdrain_07'], ['b', 'y', 'r', 'k', 'g', 'm', 'navy', 'blueviolet'], '../Plots/chip15/', 'Fresh_vs_MLC-1-7_VG1p8_VD2p0_IDS-VGS_VaS-VbD_', range(0, 128), 'Fresh vs MLC-1-7 (VG=1.8, VD=2.0)\nMLC-{1, 3, 9, 27, 81, 243, 729}ms WL pulses, IDSAT threshold = {95, 80, 65, 50, 35, 20, 5}uA, forward' , 165, ['fresh, level[0]', 'level[1]', 'level[2]', 'level[3]', 'level[4]', 'level[5]', 'level[6]', 'level[7]']) 
    #IDS_VGS(15, 33, 16, 2, 'ULVT', 128, ['../Data/chip15/Fresh_Chip15_Col33_Ids_Vgs_VAdrain_VBsource', '../Data/chip15/MLC_Chip15_Col33_1msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip15/MLC_Chip15_Col33_3msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip15/MLC_Chip15_Col33_9msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_03', '../Data/chip15/MLC_Chip15_Col33_27msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_04', '../Data/chip15/MLC_Chip15_Col33_81msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_05', '../Data/chip15/MLC_Chip15_Col33_243msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_06', '../Data/chip15/MLC_Chip15_Col33_729msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_07'], ['b', 'y', 'r', 'k', 'g', 'm', 'navy', 'blueviolet'], '../Plots/chip15/', 'Fresh_vs_MLC-1-7_VG1p8_VD2p0_IDS-VGS_VaD-VbS_', range(0, 128), 'Fresh vs MLC-1-7 (VG=1.8, VD=2.0)\nMLC-{1, 3, 9, 27, 81, 243, 729}ms WL pulses, IDSAT threshold = {95, 80, 65, 50, 35, 20, 5}uA, reversed', 165, ['fresh, level[0]', 'level[1]', 'level[2]', 'level[3]', 'level[4]', 'level[5]', 'level[6]', 'level[7]']) 


    #hist_IDS_VGS(0, 15, 33, 16, 2, 'ULVT', 128, ['../Data/chip15/Fresh_Chip15_Col33_Ids_Vgs_VAdrain_VBsource', '../Data/chip15/MLC_Chip15_Col33_1msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip15/MLC_Chip15_Col33_3msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip15/MLC_Chip15_Col33_9msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_03', '../Data/chip15/MLC_Chip15_Col33_27msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_04', '../Data/chip15/MLC_Chip15_Col33_81msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_05', '../Data/chip15/MLC_Chip15_Col33_243msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_06', '../Data/chip15/MLC_Chip15_Col33_729msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_07'], ['b', 'y', 'r', 'k', 'g', 'm', 'navy', 'blueviolet'], '../Plots/chip15/', 'Hist-IDSAT_MLC-1-7_reverse-read_', range(0, 128), 'MLC programming {1, 3, 9, 27, 81, 243, 729}ms pulses, VGS=1.8, VDS=2.0 for level[1-7]\nhistogram of read-IDSAT (VGS=VDS=0.8V)', 0, 165, 0, 165, 1000)
    #
    #t_label = []
    #for t in np.arange(0, 0.002*(71) + 0.0001, 0.002):
    #    t_label.append(str(t))
    #
    ##MLC_IDSAT_algorithm_rv1(14, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, [21], [71], [0.002], 1, np.arange(0, 0.002*(71)+0.0001, 0.002), t_label, ['../Data/chip14/MLC_programming_Chip14_Col33_2msPULSE_VG1p8_VD2p0_VAsource_VBdrain_01'], '../Plots/chip14/', 'VG1p8_VD2p0', '_rv1_cycle01_row-21', Imin=82, Imax=142)

    #for row_start in np.arange(0, 128):
    #    MLC_IDSAT_algorithm_rv1(14, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, [row_start], [71], [0.002], 1, np.arange(0, 0.002*(71)+0.0001, 0.002), t_label, ['../Data/chip14/MLC_programming_Chip14_Col33_2msPULSE_VG1p8_VD2p0_VAsource_VBdrain_01'], '../Plots/chip14/', 'VG1p8_VD2p0', '_rv1_cycle01_row_'+str(row_start).zfill(3), Imin=80, Imax=142)

    #MLC_IDSAT_algorithm_rv1(14, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, range(0, 128), [71], [0.002], 1, np.arange(0, 0.002*(71)+0.0001, 0.002), t_label, ['../Data/chip14/MLC_programming_Chip14_Col33_2msPULSE_VG1p8_VD2p0_VAsource_VBdrain_01'], '../Plots/chip14/', 'VG1p8_VD2p0', '_rv1_cycle01', Imin=80, Imax=142)

    #MLC_IDSAT_characterization(11, 30, 16, 2, 'ULVT', 1.8, 1.7, 32, range(0, 16) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD1p7', '_cycle01', 50, 160, 1)
    #MLC_IDSAT_characterization(11, 30, 16, 2, 'ULVT', 1.8, 2.0, 32, range(16, 32), [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD2p0', '_cycle01', 20, 160, 1)

    #MLC_IDSAT_characterization(11, 30, 16, 2, 'ULVT', 1.8, 1.7, 32, range(0, 16) , [40, 80], [0.01, 0.01], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col30_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p8_VD1p7', '_cycle0102', 50, 160, 1)
    #MLC_IDSAT_characterization(11, 30, 16, 2, 'ULVT', 1.8, 2.0, 32, range(16, 32), [40, 80], [0.01, 0.01], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col30_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p8_VD2p0', '_cycle0102', 20, 160, 1)

    #MLC_IDSAT_characterization(11, 30, 16, 2, 'ULVT', 1.8, 1.7, 32, range(0, 16) , [40, 80, 240], [0.01, 0.01, 0.01], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col30_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col30_HCI_240x10ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p8_VD1p7', '_cycle010203', 50, 160, 1)
    #MLC_IDSAT_characterization(11, 30, 16, 2, 'ULVT', 1.8, 2.0, 32, range(16, 32), [40, 80, 240], [0.01, 0.01, 0.01], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col30_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col30_HCI_240x10ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p8_VD2p0', '_cycle010203', 20, 160, 1)

    #MLC_IDSAT_characterization(11, 30, 16, 2, 'ULVT', 1.8, 1.7, 32, range(0, 16) , [40, 80, 240, 180], [0.01, 0.01, 0.01, 0.04], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col30_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col30_HCI_240x10ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col30_HCI_180x40ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p8_VD1p7', '_cycle01020304', 40, 160, 1)
    #MLC_IDSAT_characterization(11, 30, 16, 2, 'ULVT', 1.8, 2.0, 32, range(16, 32), [40, 80, 240, 180], [0.01, 0.01, 0.01, 0.04], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col30_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col30_HCI_240x10ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col30_HCI_180x40ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p8_VD2p0', '_cycle01020304', 10, 160, 1)

    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.5, 1.7, 128, range(0, 32)  , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p5_VD1p7', '_cycle01', 50, 160, 1)
    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.8, 1.7, 128, range(32, 64) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD1p7', '_cycle01', 50, 160, 1)
    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.5, 2.0, 128, range(64, 96) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p5_VD2p0', '_cycle01', 20, 160, 1)
    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, range(96, 128), [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD2p0', '_cycle01', 20, 160, 1)

    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.5, 1.7, 128, range(0, 32)  , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col33_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p5_VD1p7', '_cycle0102', 50, 160, 1)
    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.8, 1.7, 128, range(32, 64) , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col33_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p8_VD1p7', '_cycle0102', 50, 160, 1)
    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.5, 2.0, 128, range(64, 96) , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col33_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p5_VD2p0', '_cycle0102', 20, 160, 1)
    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, range(96, 128), [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col33_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p8_VD2p0', '_cycle0102', 20, 160, 1)

    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.5, 1.7, 128, range(0, 32)  , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col33_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col33_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p5_VD1p7', '_cycle010203', 50, 160, 1)
    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.8, 1.7, 128, range(32, 64) , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col33_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col33_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p8_VD1p7', '_cycle010203', 50, 160, 1)
    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.5, 2.0, 128, range(64, 96) , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col33_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col33_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p5_VD2p0', '_cycle010203', 20, 160, 1)
    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, range(96, 128), [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col33_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col33_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p8_VD2p0', '_cycle010203', 20, 160, 1)

    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.5, 1.7, 128, range(0, 32)  , [40, 20, 12, 36], [0.01, 0.04, 0.2, 0.2], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col33_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col33_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col33_HCI_36x200ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p5_VD1p7', '_cycle01020304', 40, 160, 1)
    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.8, 1.7, 128, range(32, 64) , [40, 20, 12, 36], [0.01, 0.04, 0.2, 0.2], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col33_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col33_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col33_HCI_36x200ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p8_VD1p7', '_cycle01020304', 40, 160, 1)
    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.5, 2.0, 128, range(64, 96) , [40, 20, 12, 36], [0.01, 0.04, 0.2, 0.2], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col33_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col33_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col33_HCI_36x200ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p5_VD2p0', '_cycle01020304', 10, 160, 1)
    #MLC_IDSAT_characterization(11, 33, 16, 2, 'ULVT', 1.8, 2.0, 128, range(96, 128), [40, 20, 12, 36], [0.01, 0.04, 0.2, 0.2], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col33_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col33_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col33_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col33_HCI_36x200ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p8_VD2p0', '_cycle01020304', 10, 160, 1)

    #MLC_IDSAT_characterization(11, 18, 36, 2, 'ULVT', 1.8, 2.0, 32, range(0, 16) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD2p0', '_cycle01', 50, 125, 1)
    #MLC_IDSAT_characterization(11, 18, 36, 2, 'ULVT', 1.8, 2.4, 32, range(16, 32), [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD2p4', '_cycle01', 20, 125, 1)

    #MLC_IDSAT_characterization(11, 18, 36, 2, 'ULVT', 1.8, 2.0, 32, range(0, 16) , [40, 80], [0.01, 0.01], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col18_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p8_VD2p0', '_cycle0102', 50, 125, 1)
    #MLC_IDSAT_characterization(11, 18, 36, 2, 'ULVT', 1.8, 2.4, 32, range(16, 32), [40, 80], [0.01, 0.01], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col18_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p8_VD2p4', '_cycle0102', 20, 125, 1)

    #MLC_IDSAT_characterization(11, 18, 36, 2, 'ULVT', 1.8, 2.0, 32, range(0, 16) , [40, 80, 240], [0.01, 0.01, 0.01], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col18_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col18_HCI_240x10ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p8_VD2p0', '_cycle010203', 50, 125, 1)
    #MLC_IDSAT_characterization(11, 18, 36, 2, 'ULVT', 1.8, 2.4, 32, range(16, 32), [40, 80, 240], [0.01, 0.01, 0.01], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col18_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col18_HCI_240x10ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p8_VD2p4', '_cycle010203', 20, 125, 1)

    #MLC_IDSAT_characterization(11, 18, 36, 2, 'ULVT', 1.8, 2.0, 32, range(0, 16) , [40, 80, 240, 180], [0.01, 0.01, 0.01, 0.04], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col18_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col18_HCI_240x10ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col18_HCI_180x40ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p8_VD2p0', '_cycle01020304', 40, 125, 1)
    #MLC_IDSAT_characterization(11, 18, 36, 2, 'ULVT', 1.8, 2.4, 32, range(16, 32), [40, 80, 240, 180], [0.01, 0.01, 0.01, 0.04], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col18_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col18_HCI_240x10ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col18_HCI_180x40ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p8_VD2p4', '_cycle01020304', 10, 125, 1)

    #MLC_IDSAT_characterization(11, 24, 20, 2, 'ULVT', 1.8, 1.8, 32, range(0, 16) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD1p8', '_cycle01', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 24, 20, 2, 'ULVT', 1.8, 2.2, 32, range(16, 32), [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD2p2', '_cycle01', 15, 150, 1)

    #MLC_IDSAT_characterization(11, 24, 20, 2, 'ULVT', 1.8, 1.8, 32, range(0, 16) , [40, 80], [0.01, 0.01], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col24_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p8_VD1p8', '_cycle0102', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 24, 20, 2, 'ULVT', 1.8, 2.2, 32, range(16, 32), [40, 80], [0.01, 0.01], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col24_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p8_VD2p2', '_cycle0102', 15, 150, 1)

    #MLC_IDSAT_characterization(11, 24, 20, 2, 'ULVT', 1.8, 1.8, 32, range(0, 16) , [40, 80, 240], [0.01, 0.01, 0.01], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col24_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col24_HCI_240x10ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p8_VD1p8', '_cycle010203', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 24, 20, 2, 'ULVT', 1.8, 2.2, 32, range(16, 32), [40, 80, 240], [0.01, 0.01, 0.01], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col24_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col24_HCI_240x10ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p8_VD2p2', '_cycle010203', 15, 150, 1)

    #MLC_IDSAT_characterization(11, 24, 20, 2, 'ULVT', 1.8, 1.8, 32, range(0, 16) , [40, 80, 240, 180], [0.01, 0.01, 0.01, 0.04], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col24_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col24_HCI_240x10ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col24_HCI_180x40ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p8_VD1p8', '_cycle01020304', 40, 150, 1)
    #MLC_IDSAT_characterization(11, 24, 20, 2, 'ULVT', 1.8, 2.2, 32, range(16, 32), [40, 80, 240, 180], [0.01, 0.01, 0.01, 0.04], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col24_HCI_80x10ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col24_HCI_240x10ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col24_HCI_180x40ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p8_VD2p2', '_cycle01020304',  5, 150, 1)


    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.5, 1.8, 128, range(0, 32)  , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p5_VD1p8', '_cycle01', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.8, 1.8, 128, range(32, 64) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD1p8', '_cycle01', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.5, 2.2, 128, range(64, 96) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p5_VD2p2', '_cycle01', 15, 150, 1)
    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.8, 2.2, 128, range(96, 128), [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD2p2', '_cycle01', 15, 150, 1)

    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.5, 1.8, 128, range(0, 32)  , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col27_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p5_VD1p8', '_cycle0102', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.8, 1.8, 128, range(32, 64) , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col27_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p8_VD1p8', '_cycle0102', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.5, 2.2, 128, range(64, 96) , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col27_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p5_VD2p2', '_cycle0102', 15, 150, 1)
    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.8, 2.2, 128, range(96, 128), [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col27_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p8_VD2p2', '_cycle0102', 15, 150, 1)

    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.5, 1.8, 128, range(0, 32)  , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col27_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col27_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p5_VD1p8', '_cycle010203', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.8, 1.8, 128, range(32, 64) , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col27_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col27_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p8_VD1p8', '_cycle010203', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.5, 2.2, 128, range(64, 96) , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col27_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col27_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p5_VD2p2', '_cycle010203', 15, 150, 1)
    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.8, 2.2, 128, range(96, 128), [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col27_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col27_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p8_VD2p2', '_cycle010203', 15, 150, 1)

    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.5, 1.8, 128, range(0, 32)  , [40, 20, 12, 36], [0.01, 0.04, 0.2, 0.2], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col27_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col27_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col27_HCI_36x200ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p5_VD1p8', '_cycle01020304', 40, 150, 1)
    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.8, 1.8, 128, range(32, 64) , [40, 20, 12, 36], [0.01, 0.04, 0.2, 0.2], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col27_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col27_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col27_HCI_36x200ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p8_VD1p8', '_cycle01020304', 20, 150, 1)
    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.5, 2.2, 128, range(64, 96) , [40, 20, 12, 36], [0.01, 0.04, 0.2, 0.2], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col27_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col27_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col27_HCI_36x200ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p5_VD2p2', '_cycle01020304',  5, 150, 1)
    #MLC_IDSAT_characterization(11, 27, 20, 2, 'ULVT', 1.8, 2.2, 128, range(96, 128), [40, 20, 12, 36], [0.01, 0.04, 0.2, 0.2], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col27_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col27_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col27_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col27_HCI_36x200ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p8_VD2p2', '_cycle01020304',  5, 150, 1)

    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.5, 1.8, 128, range(0, 32)  , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p5_VD1p8', '_cycle01', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.8, 1.8, 128, range(32, 64) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD1p8', '_cycle01', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.5, 2.2, 128, range(64, 96) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p5_VD2p2', '_cycle01', 15, 150, 1)
    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.8, 2.2, 128, range(96, 128), [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD2p2', '_cycle01', 15, 150, 1)

    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.5, 1.8, 128, range(0, 32)  , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col28_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p5_VD1p8', '_cycle0102', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.8, 1.8, 128, range(32, 64) , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col28_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p8_VD1p8', '_cycle0102', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.5, 2.2, 128, range(64, 96) , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col28_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p5_VD2p2', '_cycle0102', 15, 150, 1)
    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.8, 2.2, 128, range(96, 128), [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col28_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip11/', 'VG1p8_VD2p2', '_cycle0102', 15, 150, 1)

    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.5, 1.8, 128, range(0, 32)  , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col28_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col28_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p5_VD1p8', '_cycle010203', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.8, 1.8, 128, range(32, 64) , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col28_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col28_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p8_VD1p8', '_cycle010203', 50, 150, 1)
    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.5, 2.2, 128, range(64, 96) , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col28_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col28_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p5_VD2p2', '_cycle010203', 15, 150, 1)
    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.8, 2.2, 128, range(96, 128), [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col28_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col28_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip11/', 'VG1p8_VD2p2', '_cycle010203', 15, 150, 1)

    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.5, 1.8, 128, range(0, 32)  , [40, 20, 12, 36], [0.01, 0.04, 0.2, 0.2], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col28_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col28_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col28_HCI_36x200ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p5_VD1p8', '_cycle01020304', 40, 150, 1)
    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.8, 1.8, 128, range(32, 64) , [40, 20, 12, 36], [0.01, 0.04, 0.2, 0.2], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col28_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col28_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col28_HCI_36x200ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p8_VD1p8', '_cycle01020304', 20, 150, 1)
    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.5, 2.2, 128, range(64, 96) , [40, 20, 12, 36], [0.01, 0.04, 0.2, 0.2], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col28_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col28_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col28_HCI_36x200ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p5_VD2p2', '_cycle01020304',  5, 150, 1)
    #MLC_IDSAT_characterization(11, 28, 20, 2, 'LVT', 1.8, 2.2, 128, range(96, 128), [40, 20, 12, 36], [0.01, 0.04, 0.2, 0.2], 4, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4, 4.8, 12], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6', '3.6', '10.8'], ['../Data/chip11/Chip11_Col28_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip11/Chip11_Col28_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip11/Chip11_Col28_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03', '../Data/chip11/Chip11_Col28_HCI_36x200ms_stress_VG_ConstPulse_VAsource_VBdrain_04'], '../Plots/chip11/', 'VG1p8_VD2p2', '_cycle01020304',  5, 150, 1)


    # (L, Nfin, VT_flavor, Nrow, Imax)
    col_list = [(36, 1, 'ULVT', 32 , 60 ), (36, 1, 'LVT', 32 , 50 ), (36, 1, 'SVT', 32 , 45 ),
                (36, 1, 'ULVT', 128, 60 ), (36, 1, 'LVT', 128, 50 ), (36, 1, 'SVT', 128, 45 ),
                (20, 1, 'ULVT', 32 , 75 ), (20, 1, 'LVT', 32 , 60 ), (20, 1, 'SVT', 32 , 50 ),
                (20, 1, 'ULVT', 128, 75 ), (20, 1, 'LVT', 128, 60 ), (20, 1, 'SVT', 128, 50 ),
                (16, 1, 'ULVT', 32 , 80 ), (16, 1, 'LVT', 32 , 65 ), (16, 1, 'SVT', 32 , 60 ),
                (16, 1, 'ULVT', 128, 80 ), (16, 1, 'LVT', 128, 65 ), (16, 1, 'SVT', 128, 60 ),
                (36, 2, 'ULVT', 32 , 115), (36, 2, 'LVT', 32 , 95 ), (36, 2, 'SVT', 32 , 85 ),
                (36, 2, 'ULVT', 128, 115), (36, 2, 'LVT', 128, 95 ), (36, 2, 'SVT', 128, 85 ),                
                (20, 2, 'ULVT', 32 , 135), (20, 2, 'LVT', 32 , 115), (20, 2, 'SVT', 32 , 100),
                (20, 2, 'ULVT', 128, 135), (20, 2, 'LVT', 128, 120), (20, 2, 'SVT', 128, 100),
                (16, 2, 'ULVT', 32 , 150), (16, 2, 'LVT', 32 , 125), (16, 2, 'SVT', 32 , 115),
                (16, 2, 'ULVT', 128, 150), (16, 2, 'LVT', 128, 125), (16, 2, 'SVT', 128, 115)]

    #MLC_IDSAT_algorithm_rv1(11, 21, 36, 2, 'ULVT', 1.8, 2.4, 128, range(0, 128), [135+20], [0.2], 1, np.arange(0, 0.01*16+0.0001, 0.01), '', ['../Data/chip11/MLC_programming_Chip11_Col21_2msPULSE_VG1p8_VD2p4_VAsource_VBdrain_01'], '../Plots/chip11/', 'VG1p8_VD2p4', '_rv1_cycle01_EfficientPython')

    #MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', '0.9-1.2-1.5-1.8', 2.4, 128, range(0, 128), [59+16, 72+40, 80+31, 68+23], [0.2, 0.2, 0.2, 0.2], 4, [0, 15, 15.1, 37.5, 37.6, 59.8, 59.9, 78.1], ['0', '15', '', '37.4', '', '59.6', '', '77.8'], ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG0p9_VD2p4_VAsource_VBdrain_01', '../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG1p2_VD2p4_VAsource_VBdrain_02', '../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG1p5_VD2p4_VAsource_VBdrain_03', '../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG1p8_VD2p4_VAsource_VBdrain_04'], '../Plots/chip12/', 'VG-0p9-1p2-1p5-1p8_VD2p4', '_rv1_cycle01020304')

    t_ratio_lst = [(0, 0.17), (0.16, 0.34), (0.33, 0.505), (0.495, 0.67), (0.66, 0.84), (0.83, 1)]

    #t_label = []
    #for t in np.arange(0, 0.2*(59+16) + 0.0001, 0.2):
    #    t_label.append(str(t))
    #MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', 0.9, 2.4, 128, range(0, 128), [59+16], [0.2], 1, np.arange(0, 0.2*(59+16)+0.0001, 0.2), t_label, ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG0p9_VD2p4_VAsource_VBdrain_01'], '../Plots/chip12/', 'VG0p9_VD2p4', '_rv1_cycle01')
    #for row_start in np.arange(0, 128, 8):
    #    MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', 0.9, 2.4, 128, range(row_start, row_start+8), [59+16], [0.2], 1, np.arange(0, 0.2*(59+16)+0.0001, 0.2), t_label, ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG0p9_VD2p4_VAsource_VBdrain_01'], '../Plots/chip12/', 'VG0p9_VD2p4', '_rv1_cycle01_row'+str(row_start)+'_to_'+str(row_start+7))
    #    segment=0
    #    for t_ratio in t_ratio_lst:
    #        MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', 0.9, 2.4, 128, range(row_start, row_start+8), [59+16], [0.2], 1, np.arange(0, 0.2*(59+16)+0.0001, 0.2), t_label, ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG0p9_VD2p4_VAsource_VBdrain_01'], '../Plots/chip12/', 'VG0p9_VD2p4', '_rv1_cycle01_row'+str(row_start)+'_to_'+str(row_start+7)+'_'+str(segment), [t_ratio[0]*0.2*(59+16), t_ratio[1]*0.2*(59+16)])
    #        segment += 1

    #t_label = []
    #for t in np.arange(0, 0.2*(72+40) + 0.0001, 0.2):
    #    t_label.append(str(0.2*(59+16) + t))
    #MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', 1.2, 2.4, 128, range(0, 128), [72+40], [0.2], 1, np.arange(0, 0.2*(72+40)+0.0001, 0.2), t_label, ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG1p2_VD2p4_VAsource_VBdrain_02'], '../Plots/chip12/', 'VG1p2_VD2p4', '_rv1_cycle02')
    #for row_start in np.arange(0, 128, 8):
    #    MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', 1.2, 2.4, 128, range(row_start, row_start+8), [72+40], [0.2], 1, np.arange(0, 0.2*(72+40)+0.0001, 0.2), t_label, ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG1p2_VD2p4_VAsource_VBdrain_02'], '../Plots/chip12/', 'VG1p2_VD2p4', '_rv1_cycle02_row'+str(row_start)+'_to_'+str(row_start+7))
    #    segment=0
    #    for t_ratio in t_ratio_lst:
    #        MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', 1.2, 2.4, 128, range(row_start, row_start+8), [72+40], [0.2], 1, np.arange(0, 0.2*(72+40)+0.0001, 0.2), t_label, ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG1p2_VD2p4_VAsource_VBdrain_02'], '../Plots/chip12/', 'VG1p2_VD2p4', '_rv1_cycle02_row'+str(row_start)+'_to_'+str(row_start+7)+'_'+str(segment), [t_ratio[0]*0.2*(72+40), t_ratio[1]*0.2*(72+40)])
    #        segment += 1


    #t_label = []
    #for t in np.arange(0, 0.2*(80+31) + 0.0001, 0.2):
    #    t_label.append(str(0.2*(59+16) + 0.2*(72+40) + t))
    ##MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', 1.5, 2.4, 128, range(0, 128), [80+31], [0.2], 1, np.arange(0, 0.2*(80+31)+0.0001, 0.2), t_label, ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG1p5_VD2p4_VAsource_VBdrain_03'], '../Plots/chip12/', 'VG1p5_VD2p4', '_rv1_cycle03')
    #for row_start in np.arange(0, 128, 8):
    #    MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', 1.5, 2.4, 128, range(row_start, row_start+8), [80+31], [0.2], 1, np.arange(0, 0.2*(80+31)+0.0001, 0.2), t_label, ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG1p5_VD2p4_VAsource_VBdrain_03'], '../Plots/chip12/', 'VG1p5_VD2p4', '_rv1_cycle03_row'+str(row_start)+'_to_'+str(row_start+7))
    #    segment=0
    #    for t_ratio in t_ratio_lst:
    #        MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', 1.5, 2.4, 128, range(row_start, row_start+8), [80+31], [0.2], 1, np.arange(0, 0.2*(80+31)+0.0001, 0.2), t_label, ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG1p5_VD2p4_VAsource_VBdrain_03'], '../Plots/chip12/', 'VG1p5_VD2p4', '_rv1_cycle03_row'+str(row_start)+'_to_'+str(row_start+7)+'_'+str(segment), [t_ratio[0]*0.2*(80+31), t_ratio[1]*0.2*(80+31)])
    #        segment += 1


    #t_label = []
    #for t in np.arange(0, 0.2*(68+23) + 0.0001, 0.2):
    #    t_label.append(str(0.2*(59+16) + 0.2*(72+40) + 0.2*(80+31) + t))
    #MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', 1.8, 2.4, 128, range(0, 128), [68+23], [0.2], 1, np.arange(0, 0.2*(68+23)+0.0001, 0.2), t_label, ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG1p8_VD2p4_VAsource_VBdrain_04'], '../Plots/chip12/', 'VG1p8_VD2p4', '_rv1_cycle04')
    #for row_start in np.arange(0, 128, 8):
    #    MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', 1.8, 2.4, 128, range(row_start, row_start+8), [68+23], [0.2], 1, np.arange(0, 0.2*(68+23)+0.0001, 0.2), t_label, ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG1p8_VD2p4_VAsource_VBdrain_04'], '../Plots/chip12/', 'VG1p8_VD2p4', '_rv1_cycle04_row'+str(row_start)+'_to_'+str(row_start+7))
    #    segment=0
    #    for t_ratio in t_ratio_lst:
    #        MLC_IDSAT_algorithm_rv1(12, 21, 36, 2, 'ULVT', 1.8, 2.4, 128, range(row_start, row_start+8), [68+23], [0.2], 1, np.arange(0, 0.2*(68+23)+0.0001, 0.2), t_label, ['../Data/chip12/MLC_programming_Chip12_Col21_200msPULSE_VG1p8_VD2p4_VAsource_VBdrain_04'], '../Plots/chip12/', 'VG1p8_VD2p4', '_rv1_cycle04_row'+str(row_start)+'_to_'+str(row_start+7)+'_'+str(segment), [t_ratio[0]*0.2*(68+23), t_ratio[1]*0.2*(68+23)])
    #        segment += 1


    #MLC_IDSAT_characterization(10, 18, 36, 2, 'ULVT', 1.8, 2.0, 32, range(0, 16) , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip10/Chip10_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col18_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col18_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD2p0', '_cycle010203', 38, 112)
    #MLC_IDSAT_characterization(10, 18, 36, 2, 'ULVT', 1.8, 2.4, 32, range(16, 32), [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip10/Chip10_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col18_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col18_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD2p4', '_cycle010203', 16, 110)

    #MLC_IDSAT_characterization(10, 24, 20, 2, 'ULVT', 1.8, 1.8, 32, range(0, 16) , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip10/Chip10_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col24_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col24_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD1p8', '_cycle010203', 44, 133)
    #MLC_IDSAT_characterization(10, 24, 20, 2, 'ULVT', 1.8, 2.2, 32, range(16, 32), [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip10/Chip10_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col24_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col24_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD2p2', '_cycle010203', 14, 133)

    #MLC_IDSAT_characterization(10, 30, 16, 2, 'ULVT', 1.8, 1.7, 32, range(0, 16) , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip10/Chip10_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col30_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col30_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD1p7', '_cycle010203', 50, 135)
    #MLC_IDSAT_characterization(10, 30, 16, 2, 'ULVT', 1.8, 2.0, 32, range(16, 32), [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip10/Chip10_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col30_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col30_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD2p0', '_cycle010203', 20, 140)



    #MLC_IDSAT_characterization(10, 18, 36, 2, 'ULVT', 1.8, 2.0, 32, range(0, 16) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip10/Chip10_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD2p0', '_cycle01', 38, 112)
    #MLC_IDSAT_characterization(10, 18, 36, 2, 'ULVT', 1.8, 2.4, 32, range(16, 32), [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip10/Chip10_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD2p4', '_cycle01', 16, 110)

    #MLC_IDSAT_characterization(10, 24, 20, 2, 'ULVT', 1.8, 1.8, 32, range(0, 16) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip10/Chip10_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD1p8', '_cycle01', 44, 133)
    #MLC_IDSAT_characterization(10, 24, 20, 2, 'ULVT', 1.8, 2.2, 32, range(16, 32), [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip10/Chip10_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD2p2', '_cycle01', 14, 133)

    #MLC_IDSAT_characterization(10, 30, 16, 2, 'ULVT', 1.8, 1.7, 32, range(0, 16) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip10/Chip10_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD1p7', '_cycle01', 50, 135)
    #MLC_IDSAT_characterization(10, 30, 16, 2, 'ULVT', 1.8, 2.0, 32, range(16, 32), [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip10/Chip10_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD2p0', '_cycle01', 20, 140)



    #MLC_IDSAT_characterization(10, 18, 36, 2, 'ULVT', 1.8, 2.0, 32, range(0, 16) , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip10/Chip10_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col18_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD2p0', '_cycle0102', 38, 112)
    #MLC_IDSAT_characterization(10, 18, 36, 2, 'ULVT', 1.8, 2.4, 32, range(16, 32), [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip10/Chip10_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col18_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD2p4', '_cycle0102', 16, 110)
    #                                                                                                        
    #MLC_IDSAT_characterization(10, 24, 20, 2, 'ULVT', 1.8, 1.8, 32, range(0, 16) , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip10/Chip10_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col24_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD1p8', '_cycle0102', 44, 133)
    #MLC_IDSAT_characterization(10, 24, 20, 2, 'ULVT', 1.8, 2.2, 32, range(16, 32), [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip10/Chip10_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col24_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD2p2', '_cycle0102', 14, 133)
    #                                                                                                        
    #MLC_IDSAT_characterization(10, 30, 16, 2, 'ULVT', 1.8, 1.7, 32, range(0, 16) , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip10/Chip10_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col30_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD1p7', '_cycle0102', 50, 135)
    #MLC_IDSAT_characterization(10, 30, 16, 2, 'ULVT', 1.8, 2.0, 32, range(16, 32), [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip10/Chip10_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col30_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD2p0', '_cycle0102', 20, 140)


#    for col in range(36):
#        IDS_VGS(10, col, col_list[col][0], col_list[col][1], col_list[col][2], col_list[col][3], ['../Data/chip10/Fresh_Chip10_Col'+str(col).zfill(2)+'_Ids_Vgs_VAsource_VBdrain'], ['b'], '../Plots/chip10/', 'Fresh_Ids-Vgs_VaS-VbD_', range(0, col_list[col][3]), 'Fresh IDS-VGS, forward' , col_list[col][4])
#        IDS_VGS(10, col, col_list[col][0], col_list[col][1], col_list[col][2], col_list[col][3], ['../Data/chip10/Fresh_Chip10_Col'+str(col).zfill(2)+'_Ids_Vgs_VAdrain_VBsource'], ['b'], '../Plots/chip10/', 'Fresh_Ids-Vgs_VaD-VbS_', range(0, col_list[col][3]), 'Fresh IDS-VGS, reversed', col_list[col][4])

    #IDS_VGS(15, 21, 36, 2, 'ULVT', 128, ['../Data/chip15/Fresh_Chip15_Col21_Ids_Vgs_VAsource_VBdrain'], ['b'], '../Plots/chip15/', 'Fresh_IDS-VGS_VaS-VbD_', range(0, 128), 'Fresh IDS-VGS forward' , 130, ['fresh']) 
    #IDS_VGS(15, 21, 36, 2, 'ULVT', 128, ['../Data/chip15/Fresh_Chip15_Col21_Ids_Vgs_VAdrain_VBsource'], ['b'], '../Plots/chip15/', 'Fresh_IDS-VGS_VaD-VbS_', range(0, 128), 'Fresh IDS-VGS reversed', 130, ['fresh']) 

    #IDS_VGS(15, 33, 16, 2, 'ULVT', 128, ['../Data/chip15/Fresh_Chip15_Col33_Ids_Vgs_VAsource_VBdrain'], ['b'], '../Plots/chip15/', 'Fresh_IDS-VGS_VaS-VbD_', range(0, 128), 'Fresh IDS-VGS forward' , 165, ['fresh']) 
    #IDS_VGS(15, 33, 16, 2, 'ULVT', 128, ['../Data/chip15/Fresh_Chip15_Col33_Ids_Vgs_VAdrain_VBsource'], ['b'], '../Plots/chip15/', 'Fresh_IDS-VGS_VaD-VbS_', range(0, 128), 'Fresh IDS-VGS reversed', 165, ['fresh']) 

    #IDS_VGS(11, 21, 36, 2, 'ULVT', 128, ['../Data/chip11/Fresh_Chip11_Col21_Ids_Vgs_VAsource_VBdrain', '../Data/chip11/MLC_Chip11_Col21_2msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_01', '../Data/chip11/MLC_Chip11_Col21_10msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_02', '../Data/chip11/MLC_Chip11_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_03'], ['b', 'y', 'r', 'k'], '../Plots/chip11/', 'Fresh_vs_MLC010203_VG1p8_VD2p4_IDS-VGS_VaS-VbD_', range(0, 128), 'Fresh vs MLC-1-2-3 (VG=1.8, VD=2.4)\nMLC-{1, 2, 3}: {2ms, 10ms, 40ms} WL pulse, IDSAT threshold = {80, 60, 40}uA, forward' , 135, ['fresh', 'MLC-01', 'MLC-02', 'MLC-03']) 
    #IDS_VGS(11, 21, 36, 2, 'ULVT', 128, ['../Data/chip11/Fresh_Chip11_Col21_Ids_Vgs_VAdrain_VBsource', '../Data/chip11/MLC_Chip11_Col21_2msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip11/MLC_Chip11_Col21_10msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip11/MLC_Chip11_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_03'], ['b', 'y', 'r', 'k'], '../Plots/chip11/', 'Fresh_vs_MLC010203_VG1p8_VD2p4_IDS-VGS_VaD-VbS_', range(0, 128), 'Fresh vs MLC-1-2-3 (VG=1.8, VD=2.4)\nMLC-{1, 2, 3}: {2ms, 10ms, 40ms} WL pulse, IDSAT threshold = {80, 60, 40}uA, reversed', 135, ['fresh', 'MLC-01', 'MLC-02', 'MLC-03']) 

    #hist_IDS_VGS(0, 11, 21, 36, 2, 'ULVT', 128, ['../Data/chip11/Fresh_Chip11_Col21_Ids_Vgs_VAdrain_VBsource', '../Data/chip11/MLC_Chip11_Col21_2msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip11/MLC_Chip11_Col21_10msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip11/MLC_Chip11_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_03', '../Data/chip11/MLC_Chip11_Col21_200msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_04'], ['b', 'y', 'r', 'k', 'g'], '../Plots/chip11/', 'Hist-IDSAT_MLC-rv1-01020304_reverse-read_', range(0, 128), 'MLC programming (VGS=1.8, VDS=2.4), histogram of read-IDSAT (VGS=VDS=0.8V)', 0, 136, 0, 136, 1000)

    #IDS_VGS(14, 21, 36, 2, 'ULVT', 128, ['../Data/chip14/Fresh_Chip14_Col21_Ids_Vgs_VAsource_VBdrain', '../Data/chip14/MLC_Chip14_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_01', '../Data/chip14/MLC_Chip14_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_02', '../Data/chip14/MLC_Chip14_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_03', '../Data/chip14/MLC_Chip14_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_04'], ['b', 'y', 'r', 'k', 'g'], '../Plots/chip14/', 'Fresh_vs_MLC01020304_VG1p8_VD2p4_40msPULSE_IDS-VGS_VaS-VbD_', range(0, 128), 'Fresh vs MLC-1-2-3-4 (VG=1.8, VD=2.4)\nMLC-{1, 2, 3, 4}: all using 40ms WL pulses, IDSAT threshold = {80, 60, 40, 20}uA, forward' , 130, ['fresh', 'MLC-01', 'MLC-02', 'MLC-03', 'MLC-04']) 
    #IDS_VGS(14, 21, 36, 2, 'ULVT', 128, ['../Data/chip14/Fresh_Chip14_Col21_Ids_Vgs_VAdrain_VBsource', '../Data/chip14/MLC_Chip14_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip14/MLC_Chip14_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip14/MLC_Chip14_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_03', '../Data/chip14/MLC_Chip14_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_04'], ['b', 'y', 'r', 'k', 'g'], '../Plots/chip14/', 'Fresh_vs_MLC01020304_VG1p8_VD2p4_40msPULSE_IDS-VGS_VaD-VbS_', range(0, 128), 'Fresh vs MLC-1-2-3-4 (VG=1.8, VD=2.4)\nMLC-{1, 2, 3, 4}: all using 40ms WL pulses, IDSAT threshold = {80, 60, 40, 20}uA, reversed', 130, ['fresh', 'MLC-01', 'MLC-02', 'MLC-03', 'MLC-04']) 
    #
    #hist_IDS_VGS(0, 14, 21, 36, 2, 'ULVT', 128, ['../Data/chip14/Fresh_Chip14_Col21_Ids_Vgs_VAdrain_VBsource', '../Data/chip14/MLC_Chip14_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip14/MLC_Chip14_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip14/MLC_Chip14_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_03', '../Data/chip14/MLC_Chip14_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_04'], ['b', 'y', 'r', 'k', 'g'], '../Plots/chip14/', 'Hist-IDSAT_MLC-rv1-01020304_reverse-read_', range(0, 128), 'MLC programming always using 40ms pulses, VGS=1.8, VDS=2.4 for level=1-2-3-4\nhistogram of read-IDSAT (VGS=VDS=0.8V)', 0, 130, 0, 130, 1000)


    #IDS_VGS(10, 21, 36, 2, 'ULVT', 128, ['../Data/chip10/Fresh_Chip10_Col21_Ids_Vgs_VAsource_VBdrain', '../Data/chip10/MLC_Chip10_Col21_10msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_01', '../Data/chip10/MLC_Chip10_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_02', '../Data/chip10/MLC_Chip10_Col21_200msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_MLC010203_VG1p8_VD2p4_IDS-VGS_VaS-VbD_', range(0, 128), 'Fresh vs 2-bit/4-level MLC programming (VG=1.8, VD=2.4)\nMLC-{1, 2, 3} use {10ms, 40ms, 200ms} WL pulses, IDSAT threshold = {80uA, 60uA, 40uA}, forward read' , col_list[21][4], ['fresh', 'MLC-1', 'MLC-2', 'MLC-3']) 
    #IDS_VGS(10, 21, 36, 2, 'ULVT', 128, ['../Data/chip10/Fresh_Chip10_Col21_Ids_Vgs_VAdrain_VBsource', '../Data/chip10/MLC_Chip10_Col21_10msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip10/MLC_Chip10_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip10/MLC_Chip10_Col21_200msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_MLC010203_VG1p8_VD2p4_IDS-VGS_VaD-VbS_', range(0, 128), 'Fresh vs 2-bit/4-level MLC programming (VG=1.8, VD=2.4)\nMLC-{1, 2, 3} use {10ms, 40ms, 200ms} WL pulses, IDSAT threshold = {80uA, 60uA, 40uA}, reversed read', col_list[21][4], ['fresh', 'MLC-1', 'MLC-2', 'MLC-3']) 

    #IDSAT_vs_row(10, 21, 36, 2, 'ULVT', 128, ['../Data/chip10/Fresh_Chip10_Col21_Ids_Vgs_VAsource_VBdrain', '../Data/chip10/MLC_Chip10_Col21_10msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_01', '../Data/chip10/MLC_Chip10_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_02', '../Data/chip10/MLC_Chip10_Col21_200msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'IDSAT_vs_row_Fresh_vs_MLC010203_VG1p8_VD2p4_VaS-VbD_', range(0, 128), 'Fresh vs 2-bit/4-level MLC programming (VG=1.8, VD=2.4)\nMLC-{1, 2, 3} use {10ms, 40ms, 200ms} WL pulses, IDSAT threshold = {80uA, 60uA, 40uA}, forward read' , col_list[21][4], ['fresh', 'MLC-1', 'MLC-2', 'MLC-3']) 
    #IDSAT_vs_row(10, 21, 36, 2, 'ULVT', 128, ['../Data/chip10/Fresh_Chip10_Col21_Ids_Vgs_VAdrain_VBsource', '../Data/chip10/MLC_Chip10_Col21_10msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip10/MLC_Chip10_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip10/MLC_Chip10_Col21_200msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'IDSAT_vs_row_Fresh_vs_MLC010203_VG1p8_VD2p4_VaD-VbS_', range(0, 128), 'Fresh vs 2-bit/4-level MLC programming (VG=1.8, VD=2.4)\nMLC-{1, 2, 3} use {10ms, 40ms, 200ms} WL pulses, IDSAT threshold = {80uA, 60uA, 40uA}, reversed read', col_list[21][4], ['fresh', 'MLC-1', 'MLC-2', 'MLC-3']) 

    #IDS_VGS(10, 18, 36, 2, 'ULVT', 32, ['../Data/chip10/Fresh_Chip10_Col18_Ids_Vgs_VAsource_VBdrain', '../Data/chip10/Chip10_Col18_HCI_40x10ms_Ids_Vgs_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col18_HCI_20x40ms_Ids_Vgs_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col18_HCI_12x200ms_Ids_Vgs_VAsource_VBdrain_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_HCIstress010203_VG1p8_VD2p0_IDS-VGS_VaS-VbD_', range(0, 16) , 'Fresh vs Stress (VG=1.8, VD=2.0)\ncycle01: 40x10ms, cycle02: 20x40ms, cycle03: 12x200ms, forward' , col_list[18][4], ['fresh', 'cycle01', 'cycle02', 'cycle03']) 
    #IDS_VGS(10, 18, 36, 2, 'ULVT', 32, ['../Data/chip10/Fresh_Chip10_Col18_Ids_Vgs_VAdrain_VBsource', '../Data/chip10/Chip10_Col18_HCI_40x10ms_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip10/Chip10_Col18_HCI_20x40ms_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip10/Chip10_Col18_HCI_12x200ms_Ids_Vgs_VAdrain_VBsource_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_HCIstress010203_VG1p8_VD2p0_IDS-VGS_VaD-VbS_', range(0, 16) , 'Fresh vs Stress (VG=1.8, VD=2.0)\ncycle01: 40x10ms, cycle02: 20x40ms, cycle03: 12x200ms, reversed', col_list[18][4], ['fresh', 'cycle01', 'cycle02', 'cycle03']) 

    #IDS_VGS(10, 18, 36, 2, 'ULVT', 32, ['../Data/chip10/Fresh_Chip10_Col18_Ids_Vgs_VAsource_VBdrain', '../Data/chip10/Chip10_Col18_HCI_40x10ms_Ids_Vgs_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col18_HCI_20x40ms_Ids_Vgs_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col18_HCI_12x200ms_Ids_Vgs_VAsource_VBdrain_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_HCIstress010203_VG1p8_VD2p4_IDS-VGS_VaS-VbD_', range(16, 32), 'Fresh vs Stress (VG=1.8, VD=2.4)\ncycle01: 40x10ms, cycle02: 20x40ms, cycle03: 12x200ms, forward' , col_list[18][4], ['fresh', 'cycle01', 'cycle02', 'cycle03']) 
    #IDS_VGS(10, 18, 36, 2, 'ULVT', 32, ['../Data/chip10/Fresh_Chip10_Col18_Ids_Vgs_VAdrain_VBsource', '../Data/chip10/Chip10_Col18_HCI_40x10ms_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip10/Chip10_Col18_HCI_20x40ms_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip10/Chip10_Col18_HCI_12x200ms_Ids_Vgs_VAdrain_VBsource_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_HCIstress010203_VG1p8_VD2p4_IDS-VGS_VaD-VbS_', range(16, 32), 'Fresh vs Stress (VG=1.8, VD=2.4)\ncycle01: 40x10ms, cycle02: 20x40ms, cycle03: 12x200ms, reversed', col_list[18][4], ['fresh', 'cycle01', 'cycle02', 'cycle03']) 


    #IDS_VGS(10, 24, 20, 2, 'ULVT', 32, ['../Data/chip10/Fresh_Chip10_Col24_Ids_Vgs_VAsource_VBdrain', '../Data/chip10/Chip10_Col24_HCI_40x10ms_Ids_Vgs_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col24_HCI_20x40ms_Ids_Vgs_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col24_HCI_12x200ms_Ids_Vgs_VAsource_VBdrain_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_HCIstress010203_VG1p8_VD1p8_IDS-VGS_VaS-VbD_', range(0, 16) , 'Fresh vs Stress (VG=1.8, VD=1.8)\ncycle01: 40x10ms, cycle02: 20x40ms, cycle03: 12x200ms, forward' , col_list[24][4], ['fresh', 'cycle01', 'cycle02', 'cycle03']) 
    #IDS_VGS(10, 24, 20, 2, 'ULVT', 32, ['../Data/chip10/Fresh_Chip10_Col24_Ids_Vgs_VAdrain_VBsource', '../Data/chip10/Chip10_Col24_HCI_40x10ms_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip10/Chip10_Col24_HCI_20x40ms_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip10/Chip10_Col24_HCI_12x200ms_Ids_Vgs_VAdrain_VBsource_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_HCIstress010203_VG1p8_VD1p8_IDS-VGS_VaD-VbS_', range(0, 16) , 'Fresh vs Stress (VG=1.8, VD=1.8)\ncycle01: 40x10ms, cycle02: 20x40ms, cycle03: 12x200ms, reversed', col_list[24][4], ['fresh', 'cycle01', 'cycle02', 'cycle03']) 

    #IDS_VGS(10, 24, 20, 2, 'ULVT', 32, ['../Data/chip10/Fresh_Chip10_Col24_Ids_Vgs_VAsource_VBdrain', '../Data/chip10/Chip10_Col24_HCI_40x10ms_Ids_Vgs_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col24_HCI_20x40ms_Ids_Vgs_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col24_HCI_12x200ms_Ids_Vgs_VAsource_VBdrain_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_HCIstress010203_VG1p8_VD2p2_IDS-VGS_VaS-VbD_', range(16, 32), 'Fresh vs Stress (VG=1.8, VD=2.2)\ncycle01: 40x10ms, cycle02: 20x40ms, cycle03: 12x200ms, forward' , col_list[24][4], ['fresh', 'cycle01', 'cycle02', 'cycle03']) 
    #IDS_VGS(10, 24, 20, 2, 'ULVT', 32, ['../Data/chip10/Fresh_Chip10_Col24_Ids_Vgs_VAdrain_VBsource', '../Data/chip10/Chip10_Col24_HCI_40x10ms_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip10/Chip10_Col24_HCI_20x40ms_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip10/Chip10_Col24_HCI_12x200ms_Ids_Vgs_VAdrain_VBsource_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_HCIstress010203_VG1p8_VD2p2_IDS-VGS_VaD-VbS_', range(16, 32), 'Fresh vs Stress (VG=1.8, VD=2.2)\ncycle01: 40x10ms, cycle02: 20x40ms, cycle03: 12x200ms, reversed', col_list[24][4], ['fresh', 'cycle01', 'cycle02', 'cycle03']) 


    #IDS_VGS(10, 30, 16, 2, 'ULVT', 32, ['../Data/chip10/Fresh_Chip10_Col30_Ids_Vgs_VAsource_VBdrain', '../Data/chip10/Chip10_Col30_HCI_40x10ms_Ids_Vgs_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col30_HCI_20x40ms_Ids_Vgs_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col30_HCI_12x200ms_Ids_Vgs_VAsource_VBdrain_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_HCIstress010203_VG1p8_VD1p7_IDS-VGS_VaS-VbD_', range(0, 16) , 'Fresh vs Stress (VG=1.8, VD=1.7)\ncycle01: 40x10ms, cycle02: 20x40ms, cycle03: 12x200ms, forward' , col_list[30][4], ['fresh', 'cycle01', 'cycle02', 'cycle03']) 
    #IDS_VGS(10, 30, 16, 2, 'ULVT', 32, ['../Data/chip10/Fresh_Chip10_Col30_Ids_Vgs_VAdrain_VBsource', '../Data/chip10/Chip10_Col30_HCI_40x10ms_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip10/Chip10_Col30_HCI_20x40ms_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip10/Chip10_Col30_HCI_12x200ms_Ids_Vgs_VAdrain_VBsource_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_HCIstress010203_VG1p8_VD1p7_IDS-VGS_VaD-VbS_', range(0, 16) , 'Fresh vs Stress (VG=1.8, VD=1.7)\ncycle01: 40x10ms, cycle02: 20x40ms, cycle03: 12x200ms, reversed', col_list[30][4], ['fresh', 'cycle01', 'cycle02', 'cycle03']) 

    #IDS_VGS(10, 30, 16, 2, 'ULVT', 32, ['../Data/chip10/Fresh_Chip10_Col30_Ids_Vgs_VAsource_VBdrain', '../Data/chip10/Chip10_Col30_HCI_40x10ms_Ids_Vgs_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col30_HCI_20x40ms_Ids_Vgs_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col30_HCI_12x200ms_Ids_Vgs_VAsource_VBdrain_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_HCIstress010203_VG1p8_VD2p0_IDS-VGS_VaS-VbD_', range(16, 32), 'Fresh vs Stress (VG=1.8, VD=2.0)\ncycle01: 40x10ms, cycle02: 20x40ms, cycle03: 12x200ms, forward' , col_list[30][4], ['fresh', 'cycle01', 'cycle02', 'cycle03']) 
    #IDS_VGS(10, 30, 16, 2, 'ULVT', 32, ['../Data/chip10/Fresh_Chip10_Col30_Ids_Vgs_VAdrain_VBsource', '../Data/chip10/Chip10_Col30_HCI_40x10ms_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip10/Chip10_Col30_HCI_20x40ms_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip10/Chip10_Col30_HCI_12x200ms_Ids_Vgs_VAdrain_VBsource_03'], ['b', 'y', 'r', 'k'], '../Plots/chip10/', 'Fresh_vs_HCIstress010203_VG1p8_VD2p0_IDS-VGS_VaD-VbS_', range(16, 32), 'Fresh vs Stress (VG=1.8, VD=2.0)\ncycle01: 40x10ms, cycle02: 20x40ms, cycle03: 12x200ms, reversed', col_list[30][4], ['fresh', 'cycle01', 'cycle02', 'cycle03']) 


    #Charge_Pumping_compare([9, 9], ['Fresh', 'VG=1.8, VD=2.0, 1x1.2sec'], 21, 36, 2, 'ULVT', 128, ['../Data/chip09/Fresh_Chip09_Col21_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM', '../Data/chip09/HCIstress01_1x1p2s_Chip09_Col21_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_VG1p8_VD2p0_1p2s_Charge_Pumping_VSS_WL_0_VDD_WL_1p6', 'Fresh vs Stress (VG=1.8, VD=2.0, 1x1.2s), col[21], 5MHz Charge Pumping\n', -1.6, 0, 0.1, 17, 0, 1.6, 60, [2, 2], [[5000000, 1000], [5000000, 1000]], 0, [1, 1], '1kHz', [[0], [0]])

    #Charge_Pumping(8, 21, 36, 2, 'ULVT', 128, '../Data/chip09/Fresh_Chip09_Col21_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM', '../Plots/chip09/', 'Fresh_Charge_Pumping_VSS_WL_0_VDD_WL_1p6', '', -1.6, 0, 0.1, 17, 0, 1.6, 60, 2, [5000000, 1000], ['5MHz', '1kHz'], 0, [0, 1], 1, [0], 0)


if __name__ == '__main__':
  main()
    
