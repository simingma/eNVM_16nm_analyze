
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

def IDSAT_to_state(IDSAT, Imax=160.0, Istep=0.5):
    return int(round((Imax - IDSAT)/Istep))

def E_Xt_Dynkin(x, b, k):
    #return IDSAT_to_state(126.64, Imax=160.0, Istep=0.5) + b*np.log(1 + k*x)
    return IDSAT_to_state(144.685, Imax=160.0, Istep=0.5) + b*np.log(1 + k*x)

def pure_birth_generator_matrix(k, b0, N_state_total, Istep):
    #birth = exponential

    generator = np.zeros((N_state_total, N_state_total))
    generator[-1][-1] = 0
    for n in np.arange(N_state_total - 1):
        generator[n][n+1] = b0 * np.exp(-n * Istep * k)
        generator[n][n]   = -generator[n][n+1]

    return generator

def Testing_death_birth_generator_matrix(k, b0, d0, M, N_state_total, Istep):
    #death = linear, birth-death = exponential
    #Prove mathemetical implication: IDSAT vs time log() relationship!

    generator = np.zeros((N_state_total, N_state_total))
    #generator[-1][-1] = 0
    generator[0][1] = b0
    generator[0][0] = -b0
    for n in np.arange(1, N_state_total - 1):
        #generator[n][n-1] = d0 * (n * Istep) * np.exp(-n * Istep/13.0)
        generator[n][n-1] = d0 * (n * Istep) * np.exp(-n * Istep/M)
        #generator[n][n-1] = d0 * n
        #generator[n][n+1] = b0 * np.exp(-n * Istep * k)
        generator[n][n+1] = b0 * np.exp(-n * Istep * k) + generator[n][n-1]
        generator[n][n]   = -generator[n][n+1] - generator[n][n-1]
    generator[-1][-2] = d0 * ((N_state_total-1) * Istep) * np.exp(-(N_state_total-1) * Istep/M)
    #generator[-1][-2] = d0 * (N_state_total-1)
    generator[-1][-1] = -generator[-1][-2]

    return generator

def death_birth_generator_matrix(k, b0, d0, r, N_state_total, Istep):
    #death = linear, birth-death = exponential
    #Prove mathemetical implication: IDSAT vs time log() relationship!

    generator = np.zeros((N_state_total, N_state_total))
    #generator[-1][-1] = 0
    generator[0][1] = b0
    generator[0][0] = -b0
    for n in np.arange(1, N_state_total - 1):
        generator[n][n-1] = d0 * (n * Istep) ** r
        #generator[n][n-1] = d0 * n
        #generator[n][n+1] = b0 * np.exp(-n * Istep * k)
        generator[n][n+1] = b0 * np.exp(-n * Istep * k) + generator[n][n-1]
        generator[n][n]   = -generator[n][n+1] - generator[n][n-1]
    generator[-1][-2] = d0 * (N_state_total-1) ** r
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

    """ adapted from Marcov_Chain """

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
        NO_dup_time_IDSAT_stressed = []
        NO_dup_time_stress_times = []
        stress_times_min = 1e6 #initialize to a big number
        for row in np.arange(0, Nrow):
            IDSAT_stressed.append(np.zeros(int(np.sum(N_pulse[row])+PulseCycle)))
            stress_times.append(np.arange(0, pulse_length[0]*N_pulse[row][0]+0.0001, pulse_length[0]))
            NO_dup_time_IDSAT_stressed.append(np.zeros(int(np.sum(N_pulse[row])+1)))
            NO_dup_time_stress_times.append(np.arange(0, pulse_length[0]*N_pulse[row][0]+0.0001, pulse_length[0]))
            if PulseCycle > 1:
                t_accum = pulse_length[0]*N_pulse[row][0]
                for cycle in np.arange(1, PulseCycle):
                    stress_times[row] = np.append(stress_times[row], np.arange(t_accum, t_accum + pulse_length[cycle]*N_pulse[row][cycle]+0.0001, pulse_length[cycle]))
                    NO_dup_time_stress_times[row] = np.append(NO_dup_time_stress_times[row], np.arange(t_accum + pulse_length[cycle], t_accum + pulse_length[cycle]*N_pulse[row][cycle]+0.0001, pulse_length[cycle]))
                    t_accum += pulse_length[cycle]*N_pulse[row][cycle]
            if (stress_times[row][-1] < stress_times_min):
                #stress_times_min_points = np.copy(stress_times[row])
                stress_times_min = stress_times[row][-1]
                row_min = row
        print(stress_times_min)
        print(row_min)
        # for concatenating real stressing pulses

        # different cycles calculate log_likelihood seperately, since they have different T_mat due to different pulse widths 
        # total log likelihood for MLE need to accumulate across levels
        #k_inv_sweep = np.array([15.489, 14, 13, 12])
        #b0_sweep = np.array([21207.8, 30000, 40000])
        k_inv_sweep = np.array([15.489])
        b0_sweep = np.array([21207.8])
        #d0_sweep = np.array([330.0])
        d0_sweep = np.arange(1000, 8000, 1000)
        #d0_sweep = np.arange(280.0, 340.0, 10.0)
        M_sweep = np.arange(10, 18, 1)
        #r_sweep = np.array([0.4])
        #r_sweep = np.array([-1.0, -2.0, -3.0])
        #r_sweep = np.arange(0.3, 0.8, 0.1)
        #k_inv_sweep = np.arange(10.0, 16.0, 2.0)
        #b0_sweep = np.arange(100000.0, 150000.0, 20000.0)
        #d0_sweep = np.arange(1.0, 8.0, 2.0)
        #r_sweep = np.arange(0.6, 1.6, 0.4)

        log_likelihood = np.zeros((len(k_inv_sweep), len(b0_sweep), len(d0_sweep), len(M_sweep), PulseCycle))
        #log_likelihood = np.zeros((len(k_inv_sweep), len(b0_sweep), len(d0_sweep), len(r_sweep), PulseCycle))

        sampled_states = []
        jump_mean = []
        jump_var = []
        jump_var_plus_mean = []
        jump_var_minus_mean = []
        # Procedure: "ALL_Initial_IDSAT_", then {[cell-by-cell: ('Before_*PULSE_IDSAT_' + 'Stress_*PULSE_IDSAT_'), then "ALL_Recovered_*PULSE_IDSAT_"] x "total pulse cycles"}
        relaxation = np.zeros(Nrow * PulseCycle)
        fig_h, ax_h = plt.subplots(nrows = 1, ncols = 1)
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
                NO_dup_time_IDSAT_stressed[row][0] = IDSAT[row][0]
                stress_idx = 0
                NO_dup_time_stress_idx = 1
                for c in np.arange(cycle):
                    stress_idx += 1 + N_pulse[row][c]
                    NO_dup_time_stress_idx += N_pulse[row][c]
                IDSAT_stressed[row][int(stress_idx): int(stress_idx + N_pulse[row][cycle] + 1)] = IDSAT[row][start_idx : int(start_idx + 3*N_pulse[row][cycle] + 1): 3]
                NO_dup_time_IDSAT_stressed[row][int(NO_dup_time_stress_idx): int(NO_dup_time_stress_idx + N_pulse[row][cycle])] = IDSAT[row][start_idx + 3: int(start_idx + 3*N_pulse[row][cycle] + 1): 3]


            #only ploting the last level distribution as though we can snapshot all the IDSAT immediately after each one passing the last threshold
            IDSAT_1st_pass_Ith = np.zeros(Nrow)
            IDSAT_final = np.zeros(Nrow) #final IDSAT after all cells have reached Ith (might be slightly different from the succeeding I_V curve IDSAT data?)
            for row in np.arange(0, Nrow):
                stress_idx = 0
                for c in np.arange(cycle):
                    stress_idx += 1 + N_pulse[row][c]
                IDSAT_1st_pass_Ith[row] = 1e6 * IDSAT_stressed[row][int(stress_idx + N_pulse[row][cycle])]
                IDSAT_final[row] = 1e6 * IDSAT[row][int(start_idx + write_time[cycle]*(1+1+1))]
            relaxation[cycle*Nrow: (cycle+1)*Nrow] = IDSAT_final - IDSAT_1st_pass_Ith
            #fig, ax = plt.subplots(nrows = 1, ncols = 1)
            if (1):
                Num_bins = 1000
                #n0, bins0, patches0 = plt.hist(IDSAT_1st_pass_Ith, bins=Num_bins, normed = 1, range=(Imin, Imax), orientation = 'vertical', color = 'y', alpha = 0.8, edgecolor = 'none')
                #mean0 = np.mean(IDSAT_1st_pass_Ith)
                #std_dev0 = np.std(IDSAT_1st_pass_Ith)
                #print(mean0, std_dev0)
                n1, bins1, patches1 = ax_h.hist(IDSAT_final, bins=Num_bins, normed = 1, range=(Imin, Imax), orientation = 'horizontal', color = 'crimson', alpha = 0.8, edgecolor = 'none')
                #n1, bins1, patches1 = plt.hist(IDSAT_final, bins=Num_bins, normed = 1, range=(Imin, Imax), orientation = 'vertical', color = 'crimson', alpha = 0.8, edgecolor = 'none')
                mean1 = np.mean(IDSAT_final)
                std_dev1 = np.std(IDSAT_final)
                print('experimantal data: mean=', mean1, ',std=', std_dev1)
                #ax.grid()
                plt.ylabel('IDSAT (uA)')
                #plt.xlabel('current (uA)')
                #plt.ylabel('number of cells')
                plt.xlabel('percentage of transistors')
                #plt.ylim(ymax=0.7)
                plt.xlim(xmax = 0.7)
            #plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', level '+str(cycle)+', '+VGVD_char+', pulse_width '+str(pulse_length[cycle])+'sec\nIDSAT distribution immediately after passing Ith='+str(Ith[cycle])+' (mean='+str(mean0)+' std='+str(std_dev0)+')\nvs after relaxation (mean='+str(mean1)+' std='+str(std_dev1)+')', fontsize=9)
            #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_level-'+str(cycle+1)+'_distribution_initial-vs-relaxed'+title+'.pdf')
            #plt.close()


            """
            fig, ax = plt.subplots(nrows = 1, ncols = 1)
            IDSAT_relax = IDSAT_final - IDSAT_1st_pass_Ith
            n, bins, patches = ax.hist(IDSAT_relax, bins=40, normed = False, orientation = 'vertical', color = 'r')
            plt.xlabel('IDSAT relaxation (uA)')
            plt.ylabel('number of cells')
            mean = np.mean(IDSAT_relax)
            std_dev = np.std(IDSAT_relax)
            #print(mean, std_dev)
            plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', level '+str(cycle)+', '+VGVD_char+', pulse_width '+str(pulse_length[cycle])+'sec\ndistribution of IDSAT shift due to relaxation (mean='+str(mean)+' std='+str(std_dev)+')', fontsize=7)
            plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_level-'+str(cycle+1)+'_distribution_relaxation'+title+'.pdf')
            plt.close()

            print('level '+str(cycle)+': std_initial = '+str(std_dev0)+', std_shift = '+str(std_dev)+', std_final = '+str(std_dev1)+', sqrt(var_initial + var_shift) = '+str(np.sqrt(std_dev0**2 + std_dev**2)))
            """


            N_state_total = int(round((Imax - Imin)/Istep) + 1)
            N_state_transient = int(round((Imax - Ith[cycle])/Istep) + 1)
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

            """
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
            

            #Count_tran = np.zeros((N_state_total, N_state_total))
            #for row in np.arange(0, Nrow):
            #    stress_idx = 0
            #    for c in np.arange(cycle):
            #        stress_idx += 1 + N_pulse[row][c]
            #    for t in np.arange(int(stress_idx), int(stress_idx + N_pulse[row][cycle])):
            #        S_now = IDSAT_to_state(1e6*IDSAT_stressed[row][t], Istep = Istep)
            #        S_next = IDSAT_to_state(1e6*IDSAT_stressed[row][t+1], Istep = Istep)
            #        Count_tran[S_now][S_next] += 1

            #for S_now in np.arange(N_state_total):
            #    if np.sum(Count_tran[S_now] != 0):
            #        sum_samp = 1.0 * np.sum(Count_tran[S_now])
            #        mean = 0
            #        var = 0
            #        for S_next in np.arange(N_state_total):
            #            mean += (1.0 * (S_next - S_now)) * Count_tran[S_now][S_next] / sum_samp
            #        for S_next in np.arange(N_state_total):
            #            var += ((1.0 * (S_next - S_now - mean)) ** 2) * Count_tran[S_now][S_next] / sum_samp
            #        sampled_states.append(S_now)
            #        jump_mean.append(mean/pulse_length[cycle])
            #        jump_var.append(var/pulse_length[cycle])
            #        jump_var_plus_mean.append((var + mean)/pulse_length[cycle])
            #        jump_var_minus_mean.append((var - mean)/pulse_length[cycle])




            #for (i, k_inv) in zip(np.arange(len(k_inv_sweep)), k_inv_sweep):
            #    for (j, b0) in zip(np.arange(len(b0_sweep)), b0_sweep):
            #        for (l, d0) in zip(np.arange(len(d0_sweep)), d0_sweep):
            #            #for (m, r) in zip(np.arange(len(r_sweep)), r_sweep):
            #            for (m, M) in zip(np.arange(len(M_sweep)), M_sweep):
            #                k = 1.0/k_inv
            #                #generator = death_birth_generator_matrix(k, b0, d0, r, N_state_total, Istep)
            #                generator = Testing_death_birth_generator_matrix(k, b0, d0, M, N_state_total, Istep)
            #                T_mat, P_mat = transition_matrix(generator, pulse_length[cycle], N_state_total, N_state_transient)
            #                #log_likelihood[i][j][l][m][cycle] = np.sum(np.multiply(Count_tran, np.log(T_mat)))
            #                for y in np.arange(N_state_transient):
            #                    for z in np.arange(N_state_total):
            #                    #for z in np.arange(y, N_state_total):
            #                        #if (Count_tran[y][z] != 0) and (T_mat[y][z] > 0):
            #                        if (Count_tran[y][z] != 0):
            #                            log_likelihood[i][j][l][m][cycle] += Count_tran[y][z] * np.log(T_mat[y][z])
            #                #print('cycle '+str(cycle+1), k_inv, b0, d0, r, log_likelihood[i][j][l][m][cycle])
            #                print('cycle '+str(cycle+1), k_inv, b0, d0, M, log_likelihood[i][j][l][m][cycle])


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

        ## aggregate relaxations across all levels to get an "accurate" relaxation Gaussian(?) statistics for convolution
        #fig, ax = plt.subplots(nrows = 1, ncols = 1)
        #n, bins, patches = ax.hist(relaxation, bins=40, normed = False, orientation = 'vertical', color = 'r')
        #plt.xlabel('IDSAT relaxation (uA)')
        #plt.ylabel('number of cells')
        #mean = np.mean(relaxation)
        #std_dev = np.std(relaxation)
        ##print(mean, std_dev)
        #plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+VGVD_char+', all levels, '+'distribution of IDSAT shift due to relaxation (mean='+str(mean)+' std='+str(std_dev)+')', fontsize=7)
        #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_distribution_relaxation'+title+'.pdf')
        #plt.close()

        #print('mean_shift = '+str(mean)', std_shift = '+str(std_dev))
        ## aggregate relaxations across all levels to get an "accurate" relaxation Gaussian(?) statistics for convolution
        

        #mean, = plt.plot(sampled_states, jump_mean, 'g', marker = '.')
        #var, = plt.plot(sampled_states, jump_var, 'r', marker = '.')
        #var_plus_mean, = plt.plot(sampled_states, jump_var_plus_mean, 'y', marker = '.')
        #var_minus_mean, = plt.plot(sampled_states, jump_var_minus_mean, 'm', marker = '.')
        ##plt.xlim(10, 80)
        ##plt.ylim(ymin = -5e3, ymax = 4e4)
        #plt.legend([mean, var, var_plus_mean, var_minus_mean], ['mean = bi-di', 'var = bi+di', 'var + mean = 2*bi', 'var - mean = 2*di'], loc = 'best', fontsize = 12)
        #plt.grid()
        #plt.show()

        #for (i, k_inv) in zip(np.arange(len(k_inv_sweep)), k_inv_sweep):
        #    for (j, b0) in zip(np.arange(len(b0_sweep)), b0_sweep):
        #        for (l, d0) in zip(np.arange(len(d0_sweep)), d0_sweep):
        #            #for (m, r) in zip(np.arange(len(r_sweep)), r_sweep):
        #            for (m, M) in zip(np.arange(len(M_sweep)), M_sweep):
        #                #print(k_inv, b0, d0, r, log_likelihood[i][j][l][m], sum(log_likelihood[i][j][l][m]))
        #                print(k_inv, b0, d0, M, log_likelihood[i][j][l][m], sum(log_likelihood[i][j][l][m]))

        """
        fit = []
        #I_t_mean = np.zeros(int(stress_times_min/0.001))
        t_inter = np.arange(0, stress_times_min - 0.0005, 0.001)
        I_t_mean = np.zeros(len(t_inter))
        for row in np.arange(Nrow):
            fit.append(interp1d(NO_dup_time_stress_times[row], 1e6*NO_dup_time_IDSAT_stressed[row], kind = 'cubic'))
            I_t_mean += fit[row](t_inter) / Nrow
        plt.plot(t_inter, I_t_mean)
        plt.grid()
        plt.figure(2)
        plt.semilogx(t_inter, I_t_mean)
        plt.xlim(xmin = 1e-4)
        plt.grid()
        print(I_t_mean[0])

        Xt_mean = np.zeros(len(I_t_mean))
        for i in np.arange(len(I_t_mean)):
            Xt_mean[i] = IDSAT_to_state(I_t_mean[i], Imax=160.0, Istep=0.5)
        popt, pcov = curve_fit(E_Xt_Dynkin, t_inter, Xt_mean, bounds = (0, np.inf))
        print(popt)
        B = popt[0]
        r = popt[1]
        print(B, r)
        Istep = 0.5
        k_inv = B*Istep
        k = 1/k_inv
        #E_X0 = IDSAT_to_state(126.64, Imax=160.0, Istep=0.5)
        E_X0 = IDSAT_to_state(144.685, Imax=160.0, Istep=0.5)
        b0 = r/(k*Istep) * np.exp(k*Istep*E_X0)
        print(k_inv, b0)
        plt.figure(3)
        fitted, = plt.plot(t_inter, E_Xt_Dynkin(t_inter, *popt), color = 'b', linewidth = 4)
        data, = plt.plot(t_inter, Xt_mean, color = 'r', linewidth = 2)
        plt.xlabel('programming time t (sec)')
        plt.ylabel('average state at time t: E[X(t)]')
        plt.legend([data, fitted], ['observation data', 'function fitting'], loc = 'best')
        #plt.grid()
        plt.figure(4)
        fitted, = plt.semilogx(t_inter, E_Xt_Dynkin(t_inter, *popt), color = 'b', linewidth = 4)
        data, = plt.semilogx(t_inter, Xt_mean, color = 'r', linewidth = 2)
        plt.legend([data, fitted], ['observation data', 'function fitting'], loc = 'best')
        plt.xlabel('programming time t (sec), log scale')
        plt.ylabel('average state at time t: E[X(t)]')
        plt.xlim(xmin = 1e-4)
        #plt.grid()
        print(Xt_mean[0], E_Xt_Dynkin(0, *popt))
        plt.show()
        """


        #mean_t = np.zeros(len(stress_times_min_points))
        #std_t = np.zeros(len(stress_times_min_points))
        #I_t_common = np.zeros((Nrow, len(stress_times_min_points)))
        #for i in np.arange(len(stress_times_min_points)):
        #    for row in np.arange(0, Nrow):
        #        #print(row, i)
        #        I_t_common[row][i] = IDSAT_stressed[row][i] 
        #    mean_t[i] = np.mean(I_t_common[:, i])
        #    std_t[i] = np.std(I_t_common[:, i])

        #plt.plot(stress_times_min_points, mean_t)
        #plt.figure(2)
        #plt.plot(stress_times_min_points, std_t)
        #plt.show()

        k_inv = 15.489
        b0 = 21207.8
        #d0 = 80000.0
        #r = -1.5
        #generator = death_birth_generator_matrix(1/k_inv, b0, d0, r, N_state_total, Istep)

        #generator = Testing_death_birth_generator_matrix(1/k_inv, b0, 4000.0, 13.0, N_state_total, Istep)
        generator = Testing_death_birth_generator_matrix(1/k_inv, b0, 6000.0, 14.0, N_state_total, Istep)

        Prior = Gaussian_function(np.arange(N_state_total), IDSAT_to_state(126.641, Imax=160.0, Istep=Istep), 7.039/Istep)
        #Prior = Gaussian_function(np.arange(N_state_total), IDSAT_to_state(144.685), 5.778678/Istep)
        #relaxation = Gaussian_function(np.arange(N_state_total), IDSAT_to_state(1.35, Imax=160.0, Istep=Istep), 0.95/Istep)
        relaxation_gaussian = Gaussian_function(np.arange(N_state_total), IDSAT_to_state(1.3526, Imax=160.0, Istep=Istep), 0.9733/Istep)
        for cycle in np.arange(PulseCycle):
            N_state_transient = int(round((Imax - Ith[cycle])/Istep) + 1)
            T_mat, P_mat = transition_matrix(generator, pulse_length[cycle], N_state_total, N_state_transient) 
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
            convolve_relax = np.convolve(Prior, relaxation_gaussian)
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
            if (1):
                #plt.plot(np.arange(Imax, Imin-1e-4, -Istep), Prior/Istep, 'b', linewidth = 1.7)
                plt.plot(convolve_relax[-N_state_total:]/Istep, np.arange(Imax, Imin-1e-4, -Istep), 'darkblue', linewidth = 2.8)
                #plt.plot(np.arange(Imax, Imin-1e-4, -Istep), convolve_relax[-N_state_total:]/Istep, 'k', linewidth = 1.7)
                #plt.xlim(20, 35)
                #plt.xlim(10, 105)
                plt.ylim(20, 100)
                plt.xlim(0, 1.2)
                #plt.xlim(20, 100)
            print(sum(Prior), sum(convolve_relax[-N_state_total:]))
        ax_h.set_aspect(aspect = 0.012)
        plt.xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6'], fontsize=17)
        plt.yticks([20, 30, 40, 50, 60, 70, 80, 90, 100], ['20', '30', '40', '50', '60', '70', '80', '90', '100'], fontsize=17)




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
