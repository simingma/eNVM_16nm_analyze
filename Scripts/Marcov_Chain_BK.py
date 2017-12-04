
import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.linalg as sp_lin
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
from MLC_IDSAT import *

def Gaussian_function(x, mu, sigma):
    f_x = np.exp(-(x - mu)**2/(2*sigma**2))
    return (f_x/np.sum(f_x))

def IDSAT_to_state(IDSAT, Imax=160.0, Istep=0.2):
    return int(round((Imax - IDSAT)/Istep))

def pure_birth_generator_matrix(k, b0, N_state_total, Istep):
    #birth = exponential

    generator = np.zeros((N_state_total, N_state_total))
    generator[-1][-1] = 0
    for n in np.arange(N_state_total - 1):
        generator[n][n+1] = b0 * np.exp(-n * Istep * k)
        generator[n][n]   = -generator[n][n+1]

    return generator

def death_birth_generator_matrix(k, b0, d0, N_state_total, Istep):
    #death = linear, birth-death = exponential

    generator = np.zeros((N_state_total, N_state_total))
    #generator[-1][-1] = 0
    generator[0][1] = b0
    generator[0][0] = -b0
    for n in np.arange(1, N_state_total - 1):
        generator[n][n-1] = d0 * n ** 1.4
        #generator[n][n-1] = d0 * n
        #generator[n][n+1] = b0 * np.exp(-n * Istep * k)
        generator[n][n+1] = b0 * np.exp(-n * Istep * k) + generator[n][n-1]
        generator[n][n]   = -generator[n][n+1] - generator[n][n-1]
    generator[-1][-2] = d0 * (N_state_total-1) ** 1.4
    #generator[-1][-2] = d0 * (N_state_total-1)
    generator[-1][-1] = -generator[-1][-2]

    return generator

def transition_matrix(generator, T_pulse, N_state_total, N_state_transient):

    T_mat = sp_lin.expm(T_pulse*generator)
    #for i in np.arange(N_state_total):
    #    print(i, T_mat[i][i])
    
    #print(np.sum(T_mat, axis=1))
    ##print(np.sum(T_mat, axis=0))

    P_mat = np.copy(T_mat)
    #P_mat = T_mat
    P_mat[N_state_transient: , :] = np.zeros((N_state_total-N_state_transient, N_state_total))
    P_mat[N_state_transient: , N_state_transient:] = np.identity(N_state_total-N_state_transient)
    #print(np.sum(P_mat, axis=1))

    return T_mat, P_mat
 

def Marcov_Chain_MLE(chip, col, L, Nfin, VT_flavor, VG, VD, Nrow, row_idx, 
        write_time, pulse_length, PulseCycle,  
        data_files, path_plot, VGVD_char, title, Imax, Ith, Imin, Istep):

    """ adapted from MLC_IDSAT_algorithm_rv1 """

    if os.path.isdir(path_plot) == False:
        os.mkdir(path_plot)

    data_points = 0
    for cycle in np.arange(PulseCycle):
        data_points += 1 + write_time[cycle]*(1+1+1)

    Total_PulseCycle = 0
    for cycle in np.arange(PulseCycle):
        Total_PulseCycle += write_time[cycle]


    #the variable counting the total number of pulse for each row during each programming level cycle.
    N_pulse = np.zeros((Nrow, PulseCycle))

    Next_Pulse = np.zeros((Nrow, Total_PulseCycle)) 
    Next_pulse_round2 = np.zeros((Nrow, Total_PulseCycle)) 
    round_idx = np.zeros((Nrow, Total_PulseCycle))

    Rows_remain = np.zeros(Total_PulseCycle)
    Rows_remain_round2 = np.zeros(Total_PulseCycle)
    # Procedure: "ALL_Initial_IDSAT_", then {[cell-by-cell: ('Before_*PULSE_IDSAT_' + 'Stress_*PULSE_IDSAT_'), then "ALL_Recovered_*PULSE_IDSAT_"] x "total pulse cycles"}
    for cycle in np.arange(0, PulseCycle):
        Pulse_all = []
        Rows_all = []
        # Only one file: I always stress in VAsource_VBdrain direction, 
        # and measurements in both directions are all in this one file
        f = open(data_files[cycle],'rU')
        Pulse_all.append(re.findall(r'round=(\d), Pulse_Cycle=\d+, Next_Pulse=(\d), Next_pulse_round2=(\d)',f.read()))
        f.seek(0)
        Rows_all.append(re.findall(r'Rows_remain=(\d+), Rows_remain_round2=(\d+)',f.read()))
        f.close()

        #NEED to OPTIMIZE this python code to improve efficiency!!!
        #print(len(Pulse_all), len(Pulse_all[0])) #(1, 19840)
        #print(len(Pulse_all[0][7]))              # 3 
        if (len(Pulse_all[0]) != Nrow * write_time[cycle]):
            print('data file error!\ngrabed '+str(len(Pulse_all[0]))+' Pulse,\nbut expecting '+str(Nrow * write_time[cycle])+' data\n')
        if (len(Rows_all[0]) != write_time[cycle]):
            print('data file error!\ngrabed '+str(len(Rows_all[0]))+' Rows_remain,\nbut expecting '+str(write_time[cycle])+' data\n')

        start_idx = 0
        pulse_idx = 0
        for c in np.arange(cycle):
            start_idx += 1+write_time[c]*(1+1+1)
            pulse_idx += write_time[c]

        rows_remain_tmp, rows_remain_round2_tmp = zip(*Rows_all[0])
        Rows_remain[pulse_idx: pulse_idx+write_time[cycle]] = np.int64(list(rows_remain_tmp))
        Rows_remain_round2[pulse_idx: pulse_idx+write_time[cycle]] = np.int64(list(rows_remain_round2_tmp))
        for row in np.arange(0, Nrow):
            #print(np.int8(Apply_Pulse_all[0][Nrow*pulse_idx + row : Nrow*(pulse_idx + write_time[cycle]) : Nrow]), Apply_Pulse_all[0][Nrow*pulse_idx + row : Nrow*(pulse_idx + write_time[cycle]) : Nrow], Nrow*pulse_idx + row , Nrow*(pulse_idx + write_time[cycle]))
            round_idx_tmp, Next_Pulse_tmp, Next_pulse_round2_tmp = zip(*Pulse_all[0][row : : Nrow])
            Next_Pulse[row][pulse_idx: pulse_idx+write_time[cycle]] = np.int8(list(Next_Pulse_tmp))
            Next_pulse_round2[row][pulse_idx: pulse_idx+write_time[cycle]] = np.int8(list(Next_pulse_round2_tmp))
            round_idx[row][pulse_idx: pulse_idx+write_time[cycle]] = np.int8(list(round_idx_tmp))

            # count the total number of applied pulses during this cycle for each row
            N_pulse[row][cycle] = np.sum(Next_Pulse[row][pulse_idx: pulse_idx+write_time[cycle]])

    for (meas_dir, direction) in [('reversed', 'VAdrain_VBsource')]:

        IDSAT = np.zeros((Nrow, data_points)) 

        # for concatenating real stressing pulses, customized for each row: list can hold variable length arrays:)
        IDSAT_stressed = []
        stress_times = []
        for row in np.arange(0, Nrow):
            IDSAT_stressed.append(np.zeros(int(np.sum(N_pulse[row])+PulseCycle)))
            stress_times.append(np.arange(0, pulse_length[0]*N_pulse[row][0]+0.0001, pulse_length[0]))
            if PulseCycle > 1:
                t_accum = pulse_length[0]*N_pulse[row][0]
                for cycle in np.arange(1, PulseCycle):
                    stress_times[row] = np.append(stress_times[row], np.arange(t_accum, t_accum + pulse_length[cycle]*N_pulse[row][cycle]+0.0001, pulse_length[cycle]))
                    t_accum += pulse_length[cycle]*N_pulse[row][cycle]
        # for concatenating real stressing pulses

        # Procedure: "ALL_Initial_IDSAT_", then {[cell-by-cell: ('Before_*PULSE_IDSAT_' + 'Stress_*PULSE_IDSAT_'), then "ALL_Recovered_*PULSE_IDSAT_"] x "total pulse cycles"}
        for cycle in np.arange(0, PulseCycle):
            IDSAT_all = []
            # Only one file: I always stress in VAsource_VBdrain direction, 
            # and measurements in both directions are all in this one file
            f = open(data_files[cycle],'rU')
            IDSAT_all.append(re.findall(r'_IDSAT_WL\[\d+\]_'+direction+'=(-*\d*\.\d+)',f.read()))
            f.close()

            if (len(IDSAT_all[0]) != Nrow * (1+write_time[cycle]*(1+1+1))):
                print('data file error!\ngrabed '+str(len(IDSAT_all[0]))+' IDSAT,\nbut expecting '+str(Nrow * (1+write_time[cycle]*(1+1+1)))+' data\n')

            start_idx = 0
            pulse_idx = 0
            for c in np.arange(cycle):
                start_idx += 1+write_time[c]*(1+1+1)
                pulse_idx += write_time[c]

            for row in np.arange(0, Nrow):
                IDSAT[row][start_idx] = np.float64(IDSAT_all[0][row])
                for pulse in np.arange(write_time[cycle]):
                    IDSAT[row][start_idx+1+pulse*(1+1+1): start_idx+1+pulse*(1+1+1)+2] = np.float64(IDSAT_all[0][Nrow+pulse*(Nrow*(1+1+1))+row*(1+1): Nrow+pulse*(Nrow*(1+1+1))+row*(1+1)+2])
                    IDSAT[row][start_idx+1+pulse*(1+1+1)+2] = np.float64(IDSAT_all[0][Nrow+pulse*(Nrow*(1+1+1))+Nrow*(1+1)+row])
        
            for row in np.arange(0, Nrow):
                stress_idx = 0
                for c in np.arange(cycle):
                    stress_idx += 1 + N_pulse[row][c]
                IDSAT_stressed[row][int(stress_idx): int(stress_idx + N_pulse[row][cycle] + 1)] = IDSAT[row][start_idx : int(start_idx + 3*N_pulse[row][cycle] + 1): 3]

            """
            #only ploting the last level distribution as though we can snapshot all the IDSAT immediately after each one passing the last threshold
            #can only deal with one cycle (or the last cycle)
            IDSAT_1st_pass_Ith = np.zeros(Nrow)
            IDSAT_final = np.zeros(Nrow) #final IDSAT after all cells have reached Ith (might be slightly different from the succeeding I_V curve IDSAT data?)
            for row in np.arange(0, Nrow):
                IDSAT_1st_pass_Ith[row] = 1e6 * IDSAT_stressed[row][-1]
                IDSAT_final[row] = 1e6 * IDSAT[row][-1]
            fig, ax = plt.subplots(nrows = 1, ncols = 1)
            n0, bins0, patches0 = ax.hist(IDSAT_1st_pass_Ith, bins=40, normed = False, orientation = 'vertical', color = 'b', alpha = 0.4)
            mean0 = np.mean(IDSAT_1st_pass_Ith)
            std_dev0 = np.std(IDSAT_1st_pass_Ith)
            print(mean0, std_dev0)
            n1, bins1, patches1 = ax.hist(IDSAT_final, bins=40, normed = False, orientation = 'vertical', color = 'r', alpha = 0.4)
            mean1 = np.mean(IDSAT_final)
            std_dev1 = np.std(IDSAT_final)
            print(mean1, std_dev1)
            ax.grid()
            plt.xlabel('IDSAT (uA)')
            plt.ylabel('number of cells')

            fig, ax = plt.subplots(nrows = 1, ncols = 1)
            IDSAT_relax = IDSAT_final - IDSAT_1st_pass_Ith
            n, bins, patches = ax.hist(IDSAT_relax, bins=40, normed = False, orientation = 'vertical', color = 'r')
            plt.xlabel('IDSAT relaxation (uA)')
            plt.ylabel('number of cells')
            mean = np.mean(IDSAT_relax)
            std_dev = np.std(IDSAT_relax)
            print(mean, std_dev)
            plt.show()
            """


            N_state_total = int(round((Imax - Imin)/Istep) +1)
            N_state_transient = int(round((Imax - Ith)/Istep) + 1)
            #IDSAT_vs_times = []
            #for I in np.arange(160, 90, -1):
            #    IDSAT_vs_times.append([])
            #for row in np.arange(0, Nrow):
            #    for t in np.arange(len(IDSAT_stressed[row])):
            #        IDSAT_vs_times[IDSAT_to_state(1e6*IDSAT_stressed[row][t], Imax, Istep)].append(t)
            #IDSAT_vs_time_average = np.zeros(len(np.arange(150, 90, -1)))
            #for (i, I) in zip(np.arange(len(np.arange(150, 90, -1))), np.arange(150, 90, -1)):
            #    IDSAT_vs_time_average[i] = np.mean(IDSAT_vs_times[IDSAT_to_state(I, Imax, Istep)])
            #plt.plot(IDSAT_vs_time_average, np.arange(150, 90, -1))
            #plt.grid()

            IDSAT_vs_times_1st_cross = np.zeros((Nrow, len(np.arange(160, 90, -1))))
            Samp_count = np.zeros(len(np.arange(160, 90, -1)))
            for row in np.arange(0, Nrow):
                for t in np.arange(len(IDSAT_stressed[row])):
                    if (IDSAT_vs_times_1st_cross[row][IDSAT_to_state(1e6*IDSAT_stressed[row][t], Imax, Istep)] < 1e-7):
                        IDSAT_vs_times_1st_cross[row][IDSAT_to_state(1e6*IDSAT_stressed[row][t], Imax, Istep)] = t
                        Samp_count[IDSAT_to_state(1e6*IDSAT_stressed[row][t], Imax, Istep)] += 1
            IDSAT_vs_time_1st_cross_average = np.zeros(len(np.arange(150, 95, -1)))
            for (i, I) in zip(np.arange(len(np.arange(150, 95, -1))), np.arange(150, 95, -1)):
                IDSAT_vs_time_1st_cross_average[i] = np.sum(IDSAT_vs_times_1st_cross[ : , IDSAT_to_state(I, Imax, Istep)])/Samp_count[IDSAT_to_state(I, Imax, Istep)]
            plt.plot(pulse_length[cycle] * IDSAT_vs_time_1st_cross_average, np.arange(150, 95, -1), color = 'r', marker = '.')
            plt.grid()
            popt, pcov = curve_fit(log_shift, pulse_length[cycle] * IDSAT_vs_time_1st_cross_average, np.arange(150.0, 95.0, -1.0), bounds = (0, np.inf))
            print(popt)
            plt.plot(pulse_length[cycle] * IDSAT_vs_time_1st_cross_average, log_shift(pulse_length[cycle] * IDSAT_vs_time_1st_cross_average, *popt), color = 'g')

            plt.show()
            

            """
            Count_tran = np.zeros((N_state_total, N_state_total))
            for row in np.arange(0, Nrow):
                for t in np.arange(len(IDSAT_stressed[row]) - 1):
                    S_now = IDSAT_to_state(1e6*IDSAT_stressed[row][t])
                    S_next = IDSAT_to_state(1e6*IDSAT_stressed[row][t+1])
                    Count_tran[S_now][S_next] += 1

            log_likelihood = np.zeros((len(np.arange(20000, 100000, 20000)), len(np.arange(8.0, 36.0, 3)), len(np.arange(0.3, 2.8, 0.4))))
            for (i, b0) in zip(np.arange(len(np.arange(20000, 100000, 20000))), np.arange(20000, 100000, 20000)):
                for (j, k_r) in zip(np.arange(len(np.arange(8.0, 36.0, 3))), np.arange(8.0, 36.0, 3)):
                    for (l, d0) in zip(np.arange(len(np.arange(0.3, 2.8, 0.4))), np.arange(0.3, 2.8, 0.4)):
                        k = 1.0/k_r
                        generator = death_birth_generator_matrix(k, b0, d0, N_state_total, Istep)
                        T_mat, P_mat = transition_matrix(generator, pulse_length[cycle], N_state_total, N_state_transient)
                        #log_likelihood[i][j][l] = np.sum(np.multiply(Count_tran, np.log(T_mat)))
                        for y in np.arange(N_state_total):
                            for z in np.arange(N_state_total):
                            #for z in np.arange(y, N_state_total):
                                #if (Count_tran[y][z] != 0) and (T_mat[y][z] > 0):
                                if (Count_tran[y][z] != 0):
                                    log_likelihood[i][j][l] += Count_tran[y][z] * np.log(T_mat[y][z])
                        print(b0, k_r, d0, log_likelihood[i][j][l])
            """

            #log_likelihood = np.zeros((len(np.logspace(4, 6, num=10, base=10.0)), len(np.arange(8.0, 36.0))))
            #for (i, b0) in zip(np.arange(len(np.logspace(4, 6, num=10, base=10.0))), np.logspace(4, 6, num=10, base=10.0)):
            #    for (j, k_r) in zip(np.arange(len(np.arange(8.0, 36.0))), np.arange(8.0, 36.0)):
            #        k = 1.0/k_r
            #        generator = generator_matrix(k, b0, N_state_total, Istep)
            #        T_mat, P_mat = transition_matrix(generator, pulse_length[cycle], N_state_total, N_state_transient)
            #        #log_likelihood[i][j] = np.sum(np.multiply(Count_tran, np.log(T_mat)))
            #        for y in np.arange(N_state_total):
            #            #for z in np.arange(N_state_total):
            #            for z in np.arange(y, N_state_total):
            #                #if (Count_tran[y][z] != 0) and (T_mat[y][z] > 0):
            #                if (Count_tran[y][z] != 0):
            #                    log_likelihood[i][j] += Count_tran[y][z] * np.log(T_mat[y][z])
            #        print(b0, k_r, log_likelihood[i][j])

            #fig = plt.figure()
            #ax = fig.gca(projection='3d')
            #b0  = np.logspace(4, 6, num=10, base=10.0)
            #k_r = np.arange(8.0, 36.0)
            #b0, k_r = np.meshgrid(b0, k_r)
            #surf = ax.plot_surface(b0, k_r, log_likelihood)
            #plt.show()



            #p_b_tot = np.zeros((N_state_transient)) # total birth rate at each transient state
            #p_d_tot = np.zeros((N_state_transient)) # total death rate at each transient state
            #p_self  = np.zeros((N_state_transient)) # self transition rate at each transient state
            #for s in np.arange(N_state_transient):
            #    p_b_tot[s] = np.sum(S_matrix[s][s+1:])
            #    p_d_tot[s] = np.sum(S_matrix[s][:s])
            #    p_self[s]  = S_matrix[s][s]

            #print(p_b_tot)
            #print(p_d_tot)
            #print(p_self)

            #pb, = plt.plot(np.arange(Imax, Ith-0.01, -Istep), p_b_tot, marker = '.', color = 'g')
            #pd, = plt.plot(np.arange(Imax, Ith-0.01, -Istep), p_d_tot, marker = '.', color = 'r')
            #pself, = plt.plot(np.arange(Imax, Ith-0.01, -Istep), p_self, marker = '.', color = 'y')
            #plt.legend([pb, pd, pself], ['total birth', 'total death', 'self transition'], loc='upper left')
            #plt.xlabel('current state: IDSAT (uA)')
            #plt.ylabel('probability')
            #plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', level '+str(cycle)+', '+VGVD_char+', pulse_width '+str(pulse_length[cycle])+'sec\nself transition and total birth and death rate at each transient current state', fontsize=9)
            #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_level-'+str(cycle)+'_P_self_tot-d-b'+title+'.pdf')


def Marcov_Chain(chip, col, L, Nfin, VT_flavor, VG, VD, Nrow, row_idx, 
        write_time, pulse_length, PulseCycle,  
        data_files, path_plot, VGVD_char, title, Imax, Ith, Imin, Istep):

    """ adapted from MLC_IDSAT_algorithm_rv1 """

    if os.path.isdir(path_plot) == False:
        os.mkdir(path_plot)

    data_points = 0
    for cycle in np.arange(PulseCycle):
        data_points += 1 + write_time[cycle]*(1+1+1)

    Total_PulseCycle = 0
    for cycle in np.arange(PulseCycle):
        Total_PulseCycle += write_time[cycle]


    #the variable counting the total number of pulse for each row during each programming level cycle.
    N_pulse = np.zeros((Nrow, PulseCycle))

    Next_Pulse = np.zeros((Nrow, Total_PulseCycle)) 
    Next_pulse_round2 = np.zeros((Nrow, Total_PulseCycle)) 
    round_idx = np.zeros((Nrow, Total_PulseCycle))

    Rows_remain = np.zeros(Total_PulseCycle)
    Rows_remain_round2 = np.zeros(Total_PulseCycle)
    # Procedure: "ALL_Initial_IDSAT_", then {[cell-by-cell: ('Before_*PULSE_IDSAT_' + 'Stress_*PULSE_IDSAT_'), then "ALL_Recovered_*PULSE_IDSAT_"] x "total pulse cycles"}
    for cycle in np.arange(0, PulseCycle):
        Pulse_all = []
        Rows_all = []
        # Only one file: I always stress in VAsource_VBdrain direction, 
        # and measurements in both directions are all in this one file
        f = open(data_files[cycle],'rU')
        Pulse_all.append(re.findall(r'round=(\d), Pulse_Cycle=\d+, Next_Pulse=(\d), Next_pulse_round2=(\d)',f.read()))
        f.seek(0)
        Rows_all.append(re.findall(r'Rows_remain=(\d+), Rows_remain_round2=(\d+)',f.read()))
        f.close()

        #NEED to OPTIMIZE this python code to improve efficiency!!!
        #print(len(Pulse_all), len(Pulse_all[0])) #(1, 19840)
        #print(len(Pulse_all[0][7]))              # 3 
        if (len(Pulse_all[0]) != Nrow * write_time[cycle]):
            print('data file error!\ngrabed '+str(len(Pulse_all[0]))+' Pulse,\nbut expecting '+str(Nrow * write_time[cycle])+' data\n')
        if (len(Rows_all[0]) != write_time[cycle]):
            print('data file error!\ngrabed '+str(len(Rows_all[0]))+' Rows_remain,\nbut expecting '+str(write_time[cycle])+' data\n')

        start_idx = 0
        pulse_idx = 0
        for c in np.arange(cycle):
            start_idx += 1+write_time[c]*(1+1+1)
            pulse_idx += write_time[c]

        rows_remain_tmp, rows_remain_round2_tmp = zip(*Rows_all[0])
        Rows_remain[pulse_idx: pulse_idx+write_time[cycle]] = np.int64(list(rows_remain_tmp))
        Rows_remain_round2[pulse_idx: pulse_idx+write_time[cycle]] = np.int64(list(rows_remain_round2_tmp))
        for row in np.arange(0, Nrow):
            #print(np.int8(Apply_Pulse_all[0][Nrow*pulse_idx + row : Nrow*(pulse_idx + write_time[cycle]) : Nrow]), Apply_Pulse_all[0][Nrow*pulse_idx + row : Nrow*(pulse_idx + write_time[cycle]) : Nrow], Nrow*pulse_idx + row , Nrow*(pulse_idx + write_time[cycle]))
            round_idx_tmp, Next_Pulse_tmp, Next_pulse_round2_tmp = zip(*Pulse_all[0][row : : Nrow])
            Next_Pulse[row][pulse_idx: pulse_idx+write_time[cycle]] = np.int8(list(Next_Pulse_tmp))
            Next_pulse_round2[row][pulse_idx: pulse_idx+write_time[cycle]] = np.int8(list(Next_pulse_round2_tmp))
            round_idx[row][pulse_idx: pulse_idx+write_time[cycle]] = np.int8(list(round_idx_tmp))

            # count the total number of applied pulses during this cycle for each row
            N_pulse[row][cycle] = np.sum(Next_Pulse[row][pulse_idx: pulse_idx+write_time[cycle]])

    for (meas_dir, direction) in [('reversed', 'VAdrain_VBsource')]:

        IDSAT = np.zeros((Nrow, data_points)) 

        # for concatenating real stressing pulses, customized for each row: list can hold variable length arrays:)
        IDSAT_stressed = []
        stress_times = []
        for row in np.arange(0, Nrow):
            IDSAT_stressed.append(np.zeros(int(np.sum(N_pulse[row])+PulseCycle)))
            stress_times.append(np.arange(0, pulse_length[0]*N_pulse[row][0]+0.0001, pulse_length[0]))
            if PulseCycle > 1:
                t_accum = pulse_length[0]*N_pulse[row][0]
                for cycle in np.arange(1, PulseCycle):
                    stress_times[row] = np.append(stress_times[row], np.arange(t_accum, t_accum + pulse_length[cycle]*N_pulse[row][cycle]+0.0001, pulse_length[cycle]))
                    t_accum += pulse_length[cycle]*N_pulse[row][cycle]
        # for concatenating real stressing pulses

        # Procedure: "ALL_Initial_IDSAT_", then {[cell-by-cell: ('Before_*PULSE_IDSAT_' + 'Stress_*PULSE_IDSAT_'), then "ALL_Recovered_*PULSE_IDSAT_"] x "total pulse cycles"}
        for cycle in np.arange(0, PulseCycle):
            IDSAT_all = []
            # Only one file: I always stress in VAsource_VBdrain direction, 
            # and measurements in both directions are all in this one file
            f = open(data_files[cycle],'rU')
            IDSAT_all.append(re.findall(r'_IDSAT_WL\[\d+\]_'+direction+'=(-*\d*\.\d+)',f.read()))
            f.close()

            if (len(IDSAT_all[0]) != Nrow * (1+write_time[cycle]*(1+1+1))):
                print('data file error!\ngrabed '+str(len(IDSAT_all[0]))+' IDSAT,\nbut expecting '+str(Nrow * (1+write_time[cycle]*(1+1+1)))+' data\n')

            start_idx = 0
            pulse_idx = 0
            for c in np.arange(cycle):
                start_idx += 1+write_time[c]*(1+1+1)
                pulse_idx += write_time[c]

            for row in np.arange(0, Nrow):
                IDSAT[row][start_idx] = np.float64(IDSAT_all[0][row])
                for pulse in np.arange(write_time[cycle]):
                    IDSAT[row][start_idx+1+pulse*(1+1+1): start_idx+1+pulse*(1+1+1)+2] = np.float64(IDSAT_all[0][Nrow+pulse*(Nrow*(1+1+1))+row*(1+1): Nrow+pulse*(Nrow*(1+1+1))+row*(1+1)+2])
                    IDSAT[row][start_idx+1+pulse*(1+1+1)+2] = np.float64(IDSAT_all[0][Nrow+pulse*(Nrow*(1+1+1))+Nrow*(1+1)+row])
        
            for row in np.arange(0, Nrow):
                stress_idx = 0
                for c in np.arange(cycle):
                    stress_idx += 1 + N_pulse[row][c]
                IDSAT_stressed[row][int(stress_idx): int(stress_idx + N_pulse[row][cycle] + 1)] = IDSAT[row][start_idx : int(start_idx + 3*N_pulse[row][cycle] + 1): 3]

            N_state_total = int(round((Imax - Imin)/Istep) +1)
            S_matrix = np.zeros((N_state_total, N_state_total))
            N_state_transient = int(round((Imax - Ith)/Istep) + 1)
            S_transition = []
            for s in np.arange(N_state_transient):
                S_transition.append([])
            #print(len(S_transition))
            for row in np.arange(0, Nrow):
                for t in np.arange(len(IDSAT_stressed[row]) - 1):
                    S_now = IDSAT_to_state(1e6*IDSAT_stressed[row][t])
                    S_next = IDSAT_to_state(1e6*IDSAT_stressed[row][t+1])
                    #print(S_now, S_next)
                    S_transition[S_now].append(S_next)
                    S_matrix[S_now][S_next] += 1
            #for s in np.arange(len(S_transition)):
            #    print(s, len(S_transition[s]), S_transition[s])
            for s in np.arange(len(S_transition)):
                if(len(S_transition[s]) != 0):
                    S_matrix[s] = S_matrix[s]/len(S_transition[s])
            for S_now in np.arange(N_state_transient, N_state_total):
                for S_next in np.arange(N_state_transient, N_state_total):
                    if (S_now == S_next):
                        S_matrix[S_now][S_next] = 1
                    else:
                        S_matrix[S_now][S_next] = 0
            ##for s in np.arange(N_state_total):
            ##    print(S_matrix[s])
            #S_matrix_300 = np.linalg.matrix_power(S_matrix, 300)
            ##for s in np.arange(N_state_total):
            ##    print(S_matrix_300[s])
            #P_steady_state = np.mean(S_matrix_300, axis=0)
            ##print(P_steady_state)
            p_b_tot = np.zeros((N_state_transient)) # total birth rate at each transient state
            p_d_tot = np.zeros((N_state_transient)) # total death rate at each transient state
            p_self  = np.zeros((N_state_transient)) # self transition rate at each transient state
            for s in np.arange(N_state_transient):
                p_b_tot[s] = np.sum(S_matrix[s][s+1:])
                p_d_tot[s] = np.sum(S_matrix[s][:s])
                p_self[s]  = S_matrix[s][s]

            print(p_b_tot)
            print(p_d_tot)
            print(p_self)

            pb, = plt.plot(np.arange(Imax, Ith-0.01, -Istep), p_b_tot, marker = '.', color = 'g')
            pd, = plt.plot(np.arange(Imax, Ith-0.01, -Istep), p_d_tot, marker = '.', color = 'r')
            pself, = plt.plot(np.arange(Imax, Ith-0.01, -Istep), p_self, marker = '.', color = 'y')
            plt.legend([pb, pd, pself], ['total birth', 'total death', 'self transition'], loc='upper left')
            plt.xlabel('current state: IDSAT (uA)')
            plt.ylabel('probability')
            plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', level '+str(cycle)+', '+VGVD_char+', pulse_width '+str(pulse_length[cycle])+'sec\nself transition and total birth and death rate at each transient current state', fontsize=9)
            plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_level-'+str(cycle)+'_P_self_tot-d-b'+title+'.pdf')

            #S_matrix_power_0p001 = S_matrix[0 : N_state_transient, 0 : N_state_transient]
            #S_matrix_power_0p001 = sp_lin.fractional_matrix_power(S_matrix_power_0p001, 0.001)
            #for s in np.arange(N_state_transient):
            #    print(s, S_matrix_power_0p001[s])

            #plt.plot(np.arange(Imax, Imin-0.01, -Istep), P_steady_state, marker = '.')
            #plt.xlabel('steady-state IDSAT (uA)')
            #plt.ylabel('probability')
            #plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', level '+str(cycle)+', '+VGVD_char+', pulse_width '+str(pulse_length[cycle])+'sec\nsteady state probability of occupancy, P_matrix^300', fontsize=9)
            #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_level-'+str(cycle)+'_steady-state_probability_power300'+title+'.pdf')

 
#            if(meas_dir == 'reversed'):
#                print(last_IDSAT_above_threshold_round1[cycle])
#                fig, ax = plt.subplots(nrows = 1, ncols = 1)
#                n, bins, patches = ax.hist(last_IDSAT_above_threshold_round1[cycle], bins=40, normed = False)
#                plt.ylabel('number of cells')
#                plt.xlabel('last IDSAT before crossing threshold (uA)')
#                plt.title('level-'+str(cycle+1)+', last IDSAT before crossing threshold (uA)', fontsize=8)
#                plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_level'+str(cycle+1)+'_last_IDSAT_above_threshold.pdf')
#                plt.close()
 
#                print(delta_IDSAT_cross_threshold_round1[cycle])
#                fig, ax = plt.subplots(nrows = 1, ncols = 1)
#                n, bins, patches = ax.hist(delta_IDSAT_cross_threshold_round1[cycle], bins=40, normed = False)
#                plt.ylabel('number of cells')
#                plt.xlabel('delta IDSAT when crossing the threshold (uA)')
#                plt.title('level-'+str(cycle+1)+', delta IDSAT during the pulse when crossing threshold (uA)', fontsize=8)
#                plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_level'+str(cycle+1)+'_delta_IDSAT_cross_threshold.pdf')
#                plt.close()
        

if __name__ == '__main__':
  main()
