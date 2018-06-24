
import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

def log_shift_I0(x, a, r):
    return -a*np.log(1 + r*x)

def log_shift(x, b, a, r):
    return b - a*np.log(1 + r*x)

def exp_shift(x, a, b, c):
    return a * np.exp(-b*x) + c

def power_shift(x, b, a, n):
    return b - a * x**n

def multi_col_MLC_IDSAT_characterization(chip, col_col, L, Nfin, VT_flavor, VG, VD, Nrow, col_row_idx, 
        write_time_list, pulse_length_list, PulseCycle,  
        t, t_label,
        col_data_files, colors, path_plot, VGVD_char, title, Imin=0, Imax=0):

    """ modified from MLC_IDSAT_characterization to compare multiple columns """

# missing: it's more rigorous to extract OFF_leakages data and subtracted

    if os.path.isdir(path_plot) == False:
        os.mkdir(path_plot)

    #plt.figure(figsize=(12, 12))
    legend = []
    for (pulse_length, write_time, col, row_idx, data_files, color) in zip(pulse_length_list, write_time_list, col_col, col_row_idx, col_data_files, colors):
        #for (meas_dir, direction) in [('reversed', 'VAdrain_VBsource'), ('forward', 'VAsource_VBdrain')]:
        for (meas_dir, direction) in [('reversed', 'VAdrain_VBsource')]:
            data_points = 0
            for cycle in np.arange(PulseCycle):
                data_points += 1+1+write_time[cycle]+1

            I_mean = np.zeros((1 + data_points - 3*PulseCycle))
            IDSAT = np.zeros((Nrow, data_points)) 
            #initially all IDSAT, then cell-by-cell 'fresh'+stressing IDSAT, finally recovery all IDSAT
            for cycle in np.arange(0, PulseCycle):
                IDSAT_all = []
                # Only one file: I always stress in VAsource_VBdrain direction, 
                # and measurements in both directions are all in this one file
                f = open(data_files[cycle],'rU')
                IDSAT_all.append(re.findall(r'_IDSAT_WL\[\d+\]_'+direction+'=(-*\d*\.\d+)',f.read()))
                f.close()

                if (len(IDSAT_all[0]) != Nrow*(1+1+write_time[cycle]+1)):
                    print('data file error!\ngrabed '+str(len(IDSAT_all[0]))+' IDSAT,\nbut expecting '+str(Nrow*(1+1+write_time[cycle]+1))+' data\n')

                start_idx = 0
                for c in np.arange(cycle):
                    start_idx += 1+1+write_time[c]+1

                for row in np.arange(0, Nrow):
                    IDSAT[row][start_idx] = np.float64(IDSAT_all[0][row])
                    IDSAT[row][start_idx+1: (start_idx + 1+1+write_time[cycle]+1)-1] = np.float64(IDSAT_all[0][Nrow+row*(1+write_time[cycle]): Nrow+(row+1)*(1+write_time[cycle])])
                    IDSAT[row][(start_idx + 1+1+write_time[cycle]+1)-1] = np.float64(IDSAT_all[0][Nrow+Nrow*(1+write_time[cycle])+row])
            
            if (Imin < 0.01) and (Imax < 0.01):
                Imin = 1e6*np.amin(IDSAT[row_idx])
                Imax = 1e6*np.amax(IDSAT[row_idx])

            #plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+'\nIDSAT vs stress time, VGS='+str(VG)+'V, VDS='+str(VD)+'V\nIDSAT measured at VGS=0.8, VDS=0.8, '+meas_dir+' direction', fontsize=10)
            Tstress_accum = 0
            data_points = 0
            for cycle in np.arange(PulseCycle):
                if cycle == 0:
                    time_points = np.array([0])
                    I_mean[0] = 1e6 * np.mean(IDSAT[row_idx, 1])
                time_points = np.append(time_points, np.arange(Tstress_accum+pulse_length[cycle], Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.0001, pulse_length[cycle]))
                I_mean[data_points+1 - 3*cycle : data_points+1 - 3*cycle + write_time[cycle]] = 1e6*np.mean(IDSAT[row_idx, data_points+2: data_points+(1+1+write_time[cycle]+1)-1], axis = 0)

                #plt.plot(np.arange(Tstress_accum+pulse_length[cycle], Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.0001, pulse_length[cycle]), 1e6*IDSAT[row][data_points+2: data_points+(1+1+write_time[cycle]+1)-1], color=color, linestyle='solid', marker='.')
                #for row in row_idx:
                    #plt.plot(cycle*0.4+np.append(np.arange(Tstress_accum-pulse_length[cycle], Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.0001, pulse_length[cycle]), np.array([Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.1])), 1e6*IDSAT[row][data_points: data_points+(1+1+write_time[cycle]+1)], color=color, linestyle='solid', marker='.')

                Tstress_accum += write_time[cycle]*pulse_length[cycle]
                data_points += 1+1+write_time[cycle]+1
            axe, = plt.plot(time_points, I_mean-I_mean[0], color = color[1], marker = color[2], markersize=9, alpha=1.0, linestyle = 'None')
            legend.append(axe)
            """
            popt, pcov = curve_fit(log_shift_I0, time_points, I_mean-I_mean[0], bounds=(0, np.inf))
            plt.plot(time_points, log_shift_I0(time_points, *popt), linewidth=2.4, color=color[0])
            print(popt)
            print(1/popt[0], popt[0]*popt[1])
            """

    #plt.xticks(t, t_label, rotation=30, fontsize=9)
    #plt.xticks(t, t_l=abel)
    #plt.ylim(0, 165)
    plt.ylim(-105, 0)
    #plt.xlim(-pulse_length[0], t[-1]+0.1)
    """
    plt.xticks([0.5, 1.5, 2.5, 3.5], ['0.5', '1.5', '2.5', '3.5'], fontsize=17)
    plt.yticks([-100, -80, -60, -40, -20, 0], ['-100', '-80', '-60', '-40', '-20', '0'], fontsize=17)
    VGVD_char = ['$\mathrm{L=36nm, V_{DS}=2.0V}$', '$\mathrm{L=36nm, V_{DS}=2.4V}$', '$\mathrm{L=16nm, V_{DS}=1.7V}$', '$\mathrm{L=16nm, V_{DS}=2.0V}$']
    """
    #plt.grid()
    #plt.legend([HCI_fig, PBTI_fig], ['HCI', 'PBTI: VDS=0'], loc = 'best')
    plt.xlabel('time (sec)', fontsize=17)
    plt.ylabel('IDSAT (uA)', fontsize=17)
    plt.legend(legend, VGVD_char, fontsize=12)
    #plt.subplots_adjust(bottom=0.15)
    plt.grid()
    plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col-compare_'+'VGVD-compare_I-shift_'+meas_dir+title+'.pdf')
    #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col-compare_'+'VGVD-compare_I-shift_'+meas_dir+title+'_log-curve-fit.pdf')
    plt.xscale('log')
    #plt.show()
    #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col-compare_'+'VGVD-compare_I-shift_'+meas_dir+title+'_logTime_log-curve-fit.pdf')
    plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col-compare_'+'VGVD-compare_I-shift_'+meas_dir+title+'_logTime.pdf')
    plt.close()


def MLC_IDSAT_characterization(chip, col, L, Nfin, VT_flavor, VG, VD, Nrow, row_idx, 
        write_time, pulse_length, PulseCycle,  
        t, t_label,
        data_files, path_plot, VGVD_char, title, Imin=0, Imax=0, Iprog_plt = 0):

    """ Plot IDSAT shifts over constant VG stress
    apply write_time of pulse_length VG pulse consecutively, with quick IDSAT measurments inserted in between pulses
    after these, rest for >5 minutes of recovery, followed by IDSAT and full IDS-VGS measurement
    cycle this procedure for PulseCycle time.

    VD applied to Nhci cells, followed by Npbti cells VD=0 (FN-tunneling or so-called PBTI CVS)
    cycled until the end of Nrow. """

# missing: it's more rigorous to extract OFF_leakages data and subtracted

    if os.path.isdir(path_plot) == False:
        os.mkdir(path_plot)

    for (meas_dir, direction) in [('reversed', 'VAdrain_VBsource'), ('forward', 'VAsource_VBdrain')]:
        data_points = 0
        for cycle in np.arange(PulseCycle):
            data_points += 1+1+write_time[cycle]+1

        IDSAT = np.zeros((Nrow, data_points)) 
        #initially all IDSAT, then cell-by-cell 'fresh'+stressing IDSAT, finally recovery all IDSAT
        for cycle in np.arange(0, PulseCycle):
            IDSAT_all = []
            # Only one file: I always stress in VAsource_VBdrain direction, 
            # and measurements in both directions are all in this one file
            f = open(data_files[cycle],'rU')
            IDSAT_all.append(re.findall(r'_IDSAT_WL\[\d+\]_'+direction+'=(-*\d*\.\d+)',f.read()))
            f.close()

            if (len(IDSAT_all[0]) != Nrow*(1+1+write_time[cycle]+1)):
                print('data file error!\ngrabed '+str(len(IDSAT_all[0]))+' IDSAT,\nbut expecting '+str(Nrow*(1+1+write_time[cycle]+1))+' data\n')

            start_idx = 0
            for c in np.arange(cycle):
                start_idx += 1+1+write_time[c]+1

            for row in np.arange(0, Nrow):
                IDSAT[row][start_idx] = np.float64(IDSAT_all[0][row])
                IDSAT[row][start_idx+1: (start_idx + 1+1+write_time[cycle]+1)-1] = np.float64(IDSAT_all[0][Nrow+row*(1+write_time[cycle]): Nrow+(row+1)*(1+write_time[cycle])])
                IDSAT[row][(start_idx + 1+1+write_time[cycle]+1)-1] = np.float64(IDSAT_all[0][Nrow+Nrow*(1+write_time[cycle])+row])
        
        if (Imin < 0.01) and (Imax < 0.01):
            Imin = 1e6*np.amin(IDSAT[row_idx])
            Imax = 1e6*np.amax(IDSAT[row_idx])

        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+'\nIDSAT vs stress time, VGS='+str(VG)+'V, VDS='+str(VD)+'V\nIDSAT measured at VGS=0.8, VDS=0.8, '+meas_dir+' direction', fontsize=10)
        Tstress_accum = 0
        data_points = 0
        for cycle in np.arange(PulseCycle):

            for row in row_idx:
                # A spacing of twice the last pulse_length (the longest one, e. pulse_length[-1] = 0.2s => a spacing of 0.4s ) between each cycle to leave room for plotting the pre-stress disturbance and post-stress relaxation
                plt.plot(cycle*(pulse_length[-1]*2) + np.append(np.arange(Tstress_accum-pulse_length[cycle], Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.000001, pulse_length[cycle]), np.array([Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.5*pulse_length[-1]])), 1e6*IDSAT[row][data_points: data_points+(1+1+write_time[cycle]+1)], color='r', linestyle='solid', marker='.')

            Tstress_accum += write_time[cycle]*pulse_length[cycle]
            data_points += 1+1+write_time[cycle]+1

        #plt.xticks(t, t_label, rotation=30, fontsize=9)
        plt.xticks(t, t_label)
        plt.ylim(Imin, Imax)
        plt.xlim(-pulse_length[0], t[-1]+pulse_length[-1])
        plt.grid()
        #plt.legend([HCI_fig, PBTI_fig], ['HCI', 'PBTI: VDS=0'], loc = 'best')
        plt.xlabel('time (sec)')
        plt.ylabel('IDSAT (uA)')
        #plt.subplots_adjust(bottom=0.15)
        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_'+meas_dir+title+'.pdf')
        #plt.savefig(path_plot+'IDSAT_'+str(figN).zfill(3)+'_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+meas_dir+'.pdf')
        plt.close()

    if (Iprog_plt == 1):
        #plot ID and Isub during programming pulses, to analyze their evolution with tstress/stress levels and correlation with IDSAT degradation trajectory
        data_points = 0
        for cycle in np.arange(PulseCycle):
            data_points += write_time[cycle]

        ID_prog = np.zeros((Nrow, data_points)) 
        Isub = np.zeros((Nrow, data_points)) 

        for cycle in np.arange(0, PulseCycle):
            ID_all = []
            Isub_all = []
            f = open(data_files[cycle],'rU')
            ID_all.append(re.findall(r'WL\[\d+\]_ID_program=(-*\d*\.\d+)',f.read()))
            f.seek(0)
            Isub_all.append(re.findall(r'WL\[\d+\]_Isub=(-*\d*\.\d+)',f.read()))
            f.close()

            if (len(ID_all[0]) != Nrow*(write_time[cycle])):
                print('data file error!\ngrabed '+str(len(ID_all[0]))+' ID_prog,\nbut expecting '+str(Nrow*(write_time[cycle]))+' data\n')
            if (len(Isub_all[0]) != Nrow*(write_time[cycle])):
                print('data file error!\ngrabed '+str(len(Isub_all[0]))+' Isub,\nbut expecting '+str(Nrow*(write_time[cycle]))+' data\n')

            start_idx = 0
            for c in np.arange(cycle):
                start_idx += write_time[c]

            for row in np.arange(0, Nrow):
                ID_prog[row][start_idx: start_idx + write_time[cycle]] = np.float64(ID_all[0][row*write_time[cycle]: (row+1)*write_time[cycle]])
                Isub[row][start_idx: start_idx + write_time[cycle]] = np.float64(Isub_all[0][row*write_time[cycle]: (row+1)*write_time[cycle]])
        
        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+'\nID_prog (during programming pulses) vs stress time, VGS='+str(VG)+'V, VDS='+str(VD)+'V', fontsize=10)
        Tstress_accum = 0
        data_points = 0
        for cycle in np.arange(PulseCycle):

            for row in row_idx:
                plt.plot(cycle*(pulse_length[-1]*2)+np.arange(Tstress_accum+pulse_length[cycle], Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.0001, pulse_length[cycle]), 1e6*ID_prog[row][data_points: data_points+write_time[cycle]], color='r', linestyle='solid', marker='.')

            Tstress_accum += write_time[cycle]*pulse_length[cycle]
            data_points += write_time[cycle]

        #plt.xticks(t, t_label, rotation=30, fontsize=9)
        plt.xticks(t, t_label)
        #plt.ylim(Imin, Imax)
        #plt.xlim(-pulse_length[0], t[-1]+0.1)
        plt.grid()
        #plt.legend([HCI_fig, PBTI_fig], ['HCI', 'PBTI: VDS=0'], loc = 'best')
        plt.xlabel('time (sec)')
        plt.ylabel('ID_prog (uA)')
        #plt.subplots_adjust(bottom=0.15)
        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+title+'_ID_prog.pdf')
        #plt.savefig(path_plot+'IDSAT_'+str(figN).zfill(3)+'_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+meas_dir+'.pdf')
        plt.close()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+'\nID_prog (during programming pulses) vs stress time, VGS='+str(VG)+'V, VDS='+str(VD)+'V, log(Tstress)', fontsize=10)
        Tstress_accum = 0
        data_points = 0
        for cycle in np.arange(PulseCycle):

            for row in row_idx:
                plt.plot(np.arange(Tstress_accum+pulse_length[cycle], Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.0001, pulse_length[cycle]), 1e6*ID_prog[row][data_points: data_points+write_time[cycle]], color='r', linestyle='solid', marker='.')

            Tstress_accum += write_time[cycle]*pulse_length[cycle]
            data_points += write_time[cycle]

        #plt.xticks(t, t_label, rotation=30, fontsize=9)
        plt.xticks(t, t_label)
        #plt.ylim(Imin, Imax)
        #plt.xlim(-pulse_length[0], t[-1]+0.1)
        plt.grid()
        #plt.legend([HCI_fig, PBTI_fig], ['HCI', 'PBTI: VDS=0'], loc = 'best')
        plt.xlabel('log-scale time (sec)')
        plt.ylabel('ID_prog (uA)')
        ax.set_xscale('log')
        #plt.subplots_adjust(bottom=0.15)
        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+title+'_ID_prog_logTime.pdf')
        #plt.savefig(path_plot+'IDSAT_'+str(figN).zfill(3)+'_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+meas_dir+'.pdf')
        plt.close()


        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+'\nIsub (during programming pulses) vs stress time, VGS='+str(VG)+'V, VDS='+str(VD)+'V', fontsize=10)
        Tstress_accum = 0
        data_points = 0
        for cycle in np.arange(PulseCycle):

            for row in row_idx:
                plt.plot(cycle*(pulse_length[-1]*2)+np.arange(Tstress_accum+pulse_length[cycle], Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.0001, pulse_length[cycle]), 1e6*Isub[row][data_points: data_points+write_time[cycle]], color='g', linestyle='solid', marker='.')

            Tstress_accum += write_time[cycle]*pulse_length[cycle]
            data_points += write_time[cycle]

        #plt.xticks(t, t_label, rotation=30, fontsize=9)
        plt.xticks(t, t_label)
        #plt.ylim(Imin, Imax)
        #plt.xlim(-pulse_length[0], t[-1]+0.1)
        plt.grid()
        #plt.legend([HCI_fig, PBTI_fig], ['HCI', 'PBTI: VDS=0'], loc = 'best')
        plt.xlabel('time (sec)')
        plt.ylabel('Isub (uA)')
        #plt.subplots_adjust(bottom=0.15)
        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+title+'_Isub.pdf')
        #plt.savefig(path_plot+'IDSAT_'+str(figN).zfill(3)+'_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+meas_dir+'.pdf')
        plt.close()

def MLC_IDSAT_algorithm_naivete(chip, col, L, Nfin, VT_flavor, VG, VD, Nrow, row_idx, 
        write_time, pulse_length, PulseCycle,  
        t, t_label,
        data_files, path_plot, VGVD_char, title, Imin=0, Imax=0):

    """ adapted from MLC_IDSAT_characterization. Plotting all the IDSAT data from my first naive algorithm, with all the redundant measurements for thoroughness """

    if os.path.isdir(path_plot) == False:
        os.mkdir(path_plot)

    for (meas_dir, direction) in [('reversed', 'VAdrain_VBsource'), ('forward', 'VAsource_VBdrain')]:
        data_points = 0
        for cycle in np.arange(PulseCycle):
            data_points += 1 + write_time[cycle]*(1+1+1)

        IDSAT = np.zeros((Nrow, data_points)) 

        Total_PulseCycle = 0
        for cycle in np.arange(PulseCycle):
            Total_PulseCycle += write_time[cycle]

        #the variable counting the total number of pulse for each row during each programming level cycle.
        N_pulse = np.zeros((Nrow, PulseCycle))

        Apply_Pulse = np.zeros((Nrow, Total_PulseCycle)) 
        Rows_remain = np.zeros(Total_PulseCycle)
        # Procedure: "ALL_Initial_IDSAT_", then {[cell-by-cell: ('Before_*PULSE_IDSAT_' + 'Stress_*PULSE_IDSAT_'), then "ALL_Recovered_*PULSE_IDSAT_"] x "total pulse cycles"}
        for cycle in np.arange(0, PulseCycle):

            IDSAT_all = []
            Apply_Pulse_all = []
            Rows_remain_all = []
            # Only one file: I always stress in VAsource_VBdrain direction, 
            # and measurements in both directions are all in this one file
            f = open(data_files[cycle],'rU')
            IDSAT_all.append(re.findall(r'_IDSAT_WL\[\d+\]_'+direction+'=(-*\d*\.\d+)',f.read()))
            f.close()

            f = open(data_files[cycle],'rU')
            Apply_Pulse_all.append(re.findall(r'Apply_Pulse=(\d)',f.read()))
            f.close()

            f = open(data_files[cycle],'rU')
            Rows_remain_all.append(re.findall(r'Rows_remain=(\d+)',f.read()))
            f.close()

            if (len(IDSAT_all[0]) != Nrow * (1+write_time[cycle]*(1+1+1))):
                print('data file error!\ngrabed '+str(len(IDSAT_all[0]))+' IDSAT,\nbut expecting '+str(Nrow * (1+write_time[cycle]*(1+1+1)))+' data\n')
            if (len(Apply_Pulse_all[0]) != Nrow * write_time[cycle]):
                print('data file error!\ngrabed '+str(len(Apply_Pulse_all[0]))+' Apply_Pulse,\nbut expecting '+str(Nrow * write_time[cycle])+' data\n')
            if (len(Rows_remain_all[0]) != write_time[cycle]):
                print('data file error!\ngrabed '+str(len(Rows_remain_all[0]))+' Rows_remain,\nbut expecting '+str(write_time[cycle])+' data\n')


            start_idx = 0
            pulse_idx = 0
            for c in np.arange(cycle):
                start_idx += 1+write_time[c]*(1+1+1)
                pulse_idx += write_time[c]

            Rows_remain[pulse_idx: pulse_idx+write_time[cycle]] = np.int64(Rows_remain_all[0])
            for row in np.arange(0, Nrow):
                #print(np.int8(Apply_Pulse_all[0][Nrow*pulse_idx + row : Nrow*(pulse_idx + write_time[cycle]) : Nrow]), Apply_Pulse_all[0][Nrow*pulse_idx + row : Nrow*(pulse_idx + write_time[cycle]) : Nrow], Nrow*pulse_idx + row , Nrow*(pulse_idx + write_time[cycle]))
                Apply_Pulse[row][pulse_idx: pulse_idx+write_time[cycle]] = np.int8(Apply_Pulse_all[0][row : : Nrow])

                # count the total number of applied pulses during this cycle for each row
                N_pulse[row][cycle] = np.sum(Apply_Pulse[row][pulse_idx: pulse_idx+write_time[cycle]])

                IDSAT[row][start_idx] = np.float64(IDSAT_all[0][row])
                for pulse in np.arange(write_time[cycle]):
                    IDSAT[row][start_idx+1+pulse*(1+1+1): start_idx+1+pulse*(1+1+1)+2] = np.float64(IDSAT_all[0][Nrow+pulse*(Nrow*(1+1+1))+row*(1+1): Nrow+pulse*(Nrow*(1+1+1))+row*(1+1)+2])
                    IDSAT[row][start_idx+1+pulse*(1+1+1)+2] = np.float64(IDSAT_all[0][Nrow+pulse*(Nrow*(1+1+1))+Nrow*(1+1)+row])
        
        
        #if (Imin < 0.01) and (Imax < 0.01):
            Imin = 1e6*np.amin(IDSAT[row_idx])
            Imax = 1e6*np.amax(IDSAT[row_idx])

        """
        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+'\n, naive 2-bit MLC, IDSAT evolution, VGS='+str(VG)+'V, VDS='+str(VD)+'V\nIDSAT measured at VGS=0.8, VDS=0.8, '+meas_dir+' direction', fontsize=10)
        """
        Tstress_accum = 0
        data_points = 0
        pulse_idx = 0
        if(meas_dir == 'reversed'):
            last_IDSAT_above_threshold = np.zeros((PulseCycle, Nrow))
            delta_IDSAT_cross_threshold = np.zeros((PulseCycle, Nrow))
        for cycle in np.arange(PulseCycle):

            row_axe=[]
            marker_lst = ['.', 'o', 'D', '^', 'x', '+', 's', '*']

            for row in row_idx:
                if(len(row_idx) == len(marker_lst)):
                    row_marker = marker_lst[row - row_idx[0]]
                else:
                    row_marker = '.'

                """
                axe, = plt.plot([cycle*0.1+Tstress_accum], 1e6*IDSAT[row][data_points], color='b', marker = row_marker)
                """
                if(len(row_idx) == len(marker_lst)):
                    row_axe.append(axe)
                
                for pulse in np.arange(pulse_idx, pulse_idx+write_time[cycle]): 
                    if (Apply_Pulse[row][pulse] == 1):
                        seg_color = 'r'
                        if(meas_dir == 'reversed'):
                            if (pulse_idx+write_time[cycle] == pulse+1):
                                last_IDSAT_above_threshold[cycle][row] = 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)]
                                delta_IDSAT_cross_threshold[cycle][row] = 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)] - 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)+1]
                            elif (Apply_Pulse[row][pulse+1] == 0): 
                                last_IDSAT_above_threshold[cycle][row] = 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)]
                                delta_IDSAT_cross_threshold[cycle][row] = 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)] - 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)+1]
                    if (Apply_Pulse[row][pulse] == 0):
                        seg_color = 'b'
                    """
                    plt.plot(cycle*0.1+Tstress_accum+(pulse-pulse_idx)*pulse_length[cycle] + np.array([0, pulse_length[cycle]/(1.0+1.0+1.0)]), 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)-1: data_points+1+(pulse-pulse_idx)*(1+1+1)+1], color=seg_color, linestyle='solid', linewidth = 0.4, alpha = 0.5)
                    plt.plot(cycle*0.1+Tstress_accum+(pulse-pulse_idx)*pulse_length[cycle] + np.arange(pulse_length[cycle]/(1.0+1.0+1.0), 2*pulse_length[cycle]/(1.0+1.0+1.0)+0.0001, pulse_length[cycle]/(1.0+1.0+1.0)), 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1): data_points+1+(pulse-pulse_idx+1)*(1+1+1)-1], color=seg_color, linestyle='solid', marker=row_marker, linewidth = 0.4, markersize=2)
                    plt.plot(cycle*0.1+Tstress_accum+(pulse-pulse_idx)*pulse_length[cycle] + np.arange(2*pulse_length[cycle]/(1.0+1.0+1.0), pulse_length[cycle]+0.0001, pulse_length[cycle]/(1.0+1.0+1.0)), 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)+1: data_points+1+(pulse-pulse_idx+1)*(1+1+1)], color=seg_color, linestyle='solid', linewidth = 0.4, alpha = 0.5)
                    """

            Tstress_accum += write_time[cycle]*pulse_length[cycle]
            data_points += 1+write_time[cycle]*(1+1+1)
            pulse_idx += write_time[cycle]

            if(meas_dir == 'reversed'):
                print(last_IDSAT_above_threshold[cycle])
                fig, ax = plt.subplots(nrows = 1, ncols = 1)
                n, bins, patches = ax.hist(last_IDSAT_above_threshold[cycle], bins=40, normed = False)
                plt.ylabel('number of cells')
                plt.xlabel('last IDSAT before crossing threshold (uA)')
                plt.title('level-'+str(cycle+1)+', last IDSAT before crossing threshold (uA)', fontsize=8)
                plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_level'+str(cycle+1)+'_last_IDSAT_above_threshold.pdf')
                plt.close()

                print(delta_IDSAT_cross_threshold[cycle])
                fig, ax = plt.subplots(nrows = 1, ncols = 1)
                n, bins, patches = ax.hist(delta_IDSAT_cross_threshold[cycle], bins=40, normed = False)
                plt.ylabel('number of cells')
                plt.xlabel('delta IDSAT when crossing the threshold (uA)')
                plt.title('level-'+str(cycle+1)+', delta IDSAT during the pulse when crossing threshold (uA)', fontsize=8)
                plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_level'+str(cycle+1)+'_delta_IDSAT_cross_threshold.pdf')
                plt.close()

        """
        #plt.xticks(t, t_label, rotation=30, fontsize=9)
        plt.xticks(t, t_label, fontsize = 5)
        plt.ylim(Imin, Imax)
        plt.xlim(0, Tstress_accum+0.1*(PulseCycle-1))
        plt.grid()
        #plt.legend([HCI_fig, PBTI_fig], ['HCI', 'PBTI: VDS=0'], loc = 'best')
        plt.xlabel('time (sec)')
        plt.ylabel('IDSAT (uA)')
        #plt.subplots_adjust(bottom=0.15)

        if(len(row_idx) == len(marker_lst)):
            plt.legend(row_axe, [str(row_num) for row_num in row_idx], loc='best', fontsize = 7, numpoints = 1)
                
        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_'+meas_dir+title+'.pdf')
        #plt.savefig(path_plot+'IDSAT_'+str(figN).zfill(3)+'_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+meas_dir+'.pdf')
        plt.close()
        """

    """
    if(Nrow == len(row_idx)):
        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+'\nnaive 2-bit MLC, Rows_remain vs pulse, VGS='+str(VG)+'V, VDS='+str(VD)+'V', fontsize=10)
        pulse_idx = 0
        for cycle in np.arange(PulseCycle):
            plt.plot(pulse_idx+np.arange(0, write_time[cycle]+1), np.append(np.array([Nrow]), Rows_remain[pulse_idx:pulse_idx+write_time[cycle]]), marker='.', color = 'r')
            pulse_idx += write_time[cycle]

        t = np.arange(1, Total_PulseCycle+1)
        t_label = []
        for cycle in np.arange(PulseCycle):
            for pulse in np.arange(1, write_time[cycle]+1):
                t_label.append(str(pulse))

        plt.xticks(t, t_label, fontsize = 5)
        plt.ylim(0, Nrow)
        plt.xlim(0, Total_PulseCycle)
        plt.grid()
        #plt.legend([HCI_fig, PBTI_fig], ['HCI', 'PBTI: VDS=0'], loc = 'best')
        plt.xlabel('pulse cycle')
        plt.ylabel('Number of rows remain not reaching IDSAT threshold')
        #plt.subplots_adjust(bottom=0.15)
        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+title+'_Rows_remain.pdf')
        #plt.savefig(path_plot+'IDSAT_'+str(figN).zfill(3)+'_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+meas_dir+'.pdf')
        plt.close()
        """

    """
    print(N_pulse)
    print(np.corrcoef(N_pulse, rowvar=False))
    print(np.cov(N_pulse, rowvar=False))
    for i in np.arange(0, PulseCycle-1):
        for j in np.arange(i+1, PulseCycle):
            plt.scatter(N_pulse[:, i], N_pulse[:, j])
            plt.grid()
            plt.xlabel('Number of pulse, level '+str(i+1))
            plt.ylabel('Number of pulse, level '+str(j+1))
            plt.title('Number of pulse applied to each of the '+str(Nrow)+' cells, to program level '+str(j+1)+' vs level '+str(i+1))
            plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_Npulse_level_'+str(j+1)+'_vs_'+str(i+1)+'.pdf')
            plt.close()
            """


def MLC_IDSAT_algorithm_rv1(chip, col, L, Nfin, VT_flavor, VG, VD, Nrow, row_idx, 
        write_time, pulse_length, PulseCycle,  
        t, t_label,
        data_files, path_plot, VGVD_char, title, t_range=[], Imin=0, Imax=0):

    """ adapted from MLC_IDSAT_characterization. Plotting all the IDSAT data from my first naive algorithm, with all the redundant measurements for thoroughness """

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

    print(N_pulse)
    print(np.corrcoef(N_pulse, rowvar=False))
    print(np.cov(N_pulse, rowvar=False))
    for i in np.arange(0, PulseCycle-1):
        for j in np.arange(i+1, PulseCycle):
            plt.scatter(N_pulse[:, i], N_pulse[:, j])
            plt.grid()
            plt.xlabel('Number of pulse, level '+str(i+1))
            plt.ylabel('Number of pulse, level '+str(j+1))
            plt.title('Number of pulse applied to each of the '+str(Nrow)+' cells, to program level '+str(j+1)+' vs level '+str(i+1))
            plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_Npulse_level_'+str(j+1)+'_vs_'+str(i+1)+'.pdf')
            plt.close()


    """
    for (meas_dir, direction) in [('reversed', 'VAdrain_VBsource'), ('forward', 'VAsource_VBdrain')]:

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

        ##if (Imin < 0.01) and (Imax < 0.01):
        #    Imin = 1e6*np.amin(IDSAT[row_idx])
        #    Imax = 1e6*np.amax(IDSAT[row_idx])
        
        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+'\n, rv1 2-bit MLC, IDSAT evolution, VGS='+str(VG)+'V, VDS='+str(VD)+'V\nIDSAT measured at VGS=0.8, VDS=0.8, '+meas_dir+' direction', fontsize=10)

        for row in row_idx:
            #popt, pcov = curve_fit(log_shift, xdata, ydata, bounds=([0, 8, 0], [np.inf, 12, np.inf]))
            #popt, pcov = curve_fit(log_shift, xdata, ydata, bounds=(-np.inf, np.inf))
            popt, pcov = curve_fit(log_shift_I0, stress_times[row], 1e6*(IDSAT_stressed[row] -IDSAT_stressed[row][0]) , bounds=(0, np.inf))
            rmse = np.mean((1e6*(IDSAT_stressed[row]-IDSAT_stressed[row][0]) - log_shift_I0(stress_times[row], *popt))**2)**0.5
            #popt, pcov = curve_fit(log_shift, stress_times[row], 1e6*IDSAT_stressed[row], bounds=(0, np.inf))
            #rmse = np.mean((1e6*IDSAT_stressed[row] - log_shift(stress_times[row], *popt))**2)**0.5
            #plt.plot(stress_times[row], log_shift(stress_times[row], *popt) - log_shift(stress_times[row][0], *popt) , linewidth=0.4, color='g')
            #plt.plot(stress_times[row], 1e6*(IDSAT_stressed[row] - IDSAT_stressed[row][0]), marker = '.', markersize=2, linewidth=0.4, color='r')
            plot_times = stress_times[row]
            plot_times[0] = 0.0001  # use 0.1ms as the time of the first (fresh) current in order for it to show up in the log-plot (approximation error?)
            #plt.plot(plot_times, log_shift(stress_times[row], *popt), linewidth=0.4, color='g')
            #plt.plot(plot_times, 1e6*IDSAT_stressed[row], marker = '.', markersize=2, linewidth=0.4, color='r')
            plt.plot(plot_times, 1e6*IDSAT_stressed[row][0] + log_shift_I0(stress_times[row], *popt), linewidth=0.4, color='g')
            plt.plot(plot_times, 1e6*IDSAT_stressed[row], marker = '.', markersize=2, linewidth=0.4, color='r')
            #plt.plot(stress_times[row], log_shift(stress_times[row], *popt), linewidth=0.4, color='g')
            #plt.plot(stress_times[row], 1e6*IDSAT_stressed[row], marker = '.', markersize=2, linewidth=0.4, color='r')
            print(row, rmse, popt, 1e6*IDSAT_stressed[row][0], stress_times[row][-1])
            if len(row_idx) == 1:
                plt.xlabel('real stressing time (sec) of row_'+str(row)+'\nrmse='+str(rmse)+'uA, y = '+str(1e6*IDSAT_stressed[row][0])+' - '+str(popt[0])+' * log(1 + '+str(popt[1])+' * x)', fontsize=7)
                #plt.xlabel('real stressing time (sec) of row_'+str(row)+'\nrmse='+str(rmse)+'uA, y = '+str(popt[0])+' - '+str(popt[1])+' * log(1 + '+str(popt[2])+' * x)', fontsize=7)

        if (Imin > 0.01) and (Imax > 0.01):
            plt.ylim(Imin, Imax)
        plt.grid()
        plt.ylabel('IDSAT (uA)')

        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_'+meas_dir+title+'_log-fit_multi-level-stressing.pdf')
        #plt.ylabel('delta_I = I(t_stress) - I(t=0,fresh)')
        #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_'+meas_dir+title+'_delta-I_log-fit_multi-level-stressing.pdf')
        plt.xscale('log')
        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_'+meas_dir+title+'_LogTime_log-fit_multi-level-stressing.pdf')
        #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_'+meas_dir+title+'_delta-I_LogTime_log-fit_multi-level-stressing.pdf')
        plt.close()
        """

        
    """
        Tstress_accum = 0
        data_points = 0
        pulse_idx = 0
        if(meas_dir == 'reversed'):
            last_IDSAT_above_threshold_round1 = np.zeros((PulseCycle, Nrow))
            delta_IDSAT_cross_threshold_round1 = np.zeros((PulseCycle, Nrow))

        for cycle in np.arange(PulseCycle):

            row_axe=[]
            marker_lst = ['.', 'o', 'D', '^', 'x', '+', 's', '*']

            for row in row_idx:
                if(len(row_idx) == len(marker_lst)):
                    row_marker = marker_lst[row - row_idx[0]]
                else:
                    row_marker = '.'

                
                #careful: this only works for the first level cycle
                xdata = np.arange(0, N_pulse[row][cycle]*pulse_length[cycle]+0.0001, pulse_length[cycle])
                ydata = 1e6*IDSAT[row][data_points: int(data_points+3*N_pulse[row][cycle]+1): 3]

                #popt, pcov = curve_fit(log_shift, xdata, ydata, bounds=([0, 8, 0], [np.inf, 12, np.inf]))
                popt, pcov = curve_fit(log_shift, xdata, ydata, bounds=(-np.inf, np.inf))
                #popt, pcov = curve_fit(log_shift, xdata, ydata, bounds=(0, np.inf))
                rmse = np.mean((ydata - log_shift(xdata, *popt))**2)**0.5
                plt.plot(xdata, log_shift(xdata, *popt), marker = row_marker, markersize=2, linewidth=0.4, color='g')
                print(row, rmse, popt)
                plt.text(0.3, 90, 'rmse='+str(rmse)+'uA\ny = '+str(popt[0])+' - '+str(popt[1])+' * log(1 + '+str(popt[2])+' * x)', fontsize=9)

                #popt, pcov = curve_fit(exp_shift, xdata, ydata)
                #rmse = np.mean((ydata - exp_shift(xdata, *popt))**2)**0.5
                #print(row, rmse)
                #plt.plot(xdata, exp_shift(xdata, *popt), marker = row_marker, markersize=2, linewidth=0.4, color='g')

                #popt, pcov = curve_fit(power_shift, xdata, ydata, bounds=(0, [np.inf, np.inf, 1.]))
                #plt.plot(xdata, power_shift(xdata, *popt), marker = row_marker, markersize=2, linewidth=0.4, color='g')
                #rmse = np.mean((ydata - power_shift(xdata, *popt))**2)**0.5
                #print(row, rmse, popt)
                #plt.text(0.08, 105, 'rmse='+str(rmse)+'uA\ny = '+str(popt[0])+' - '+str(popt[1])+' * x ** '+str(popt[2]), fontsize=9)
                #careful: this only works for the first level cycle
                

                axe, = plt.plot([cycle*0.1+Tstress_accum], 1e6*IDSAT[row][data_points], color='b', marker = row_marker)
                
                if(len(row_idx) == len(marker_lst)):
                    row_axe.append(axe)
                
                for pulse in np.arange(pulse_idx, pulse_idx+write_time[cycle]): 
                    if ((round_idx[row][pulse] == 1) and (Next_Pulse[row][pulse] == 1)):
                        seg_color = 'r'

                        #CAUTIOUS: not carefully considering corner cases! Doesn't work for the slowest cell if round2 is skipped!!!(chip14)
                        if ((meas_dir == 'reversed') and (Next_Pulse[row][pulse+1] == 0)):
                            last_IDSAT_above_threshold_round1[cycle][row] = 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)-1]
                            delta_IDSAT_cross_threshold_round1[cycle][row] = -1e6*IDSAT[row][data_points+1+(pulse-pulse_idx+1)*(1+1+1)-1] + 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)-1] 

                    if (((round_idx[row][pulse] == 1) and (Next_Pulse[row][pulse] == 0)) or ((round_idx[row][pulse] == 2) and (Next_pulse_round2[row][pulse] == 0))):
                        seg_color = 'b'
                    if ((round_idx[row][pulse] == 2) and (Next_pulse_round2[row][pulse] == 1)):
                        seg_color = 'y'
                    
                    plt.plot(cycle*0.1+Tstress_accum+(pulse-pulse_idx)*pulse_length[cycle] + np.array([0, pulse_length[cycle]/(1.0+1.0+1.0)]), 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)-1: data_points+1+(pulse-pulse_idx)*(1+1+1)+1], color=seg_color, linestyle='solid', linewidth = 0.4, alpha = 0.5)
                    plt.plot(cycle*0.1+Tstress_accum+(pulse-pulse_idx)*pulse_length[cycle] + np.arange(pulse_length[cycle]/(1.0+1.0+1.0), 2*pulse_length[cycle]/(1.0+1.0+1.0)+0.0001, pulse_length[cycle]/(1.0+1.0+1.0)), 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1): data_points+1+(pulse-pulse_idx+1)*(1+1+1)-1], color=seg_color, linestyle='solid', marker=row_marker, linewidth = 0.4, markersize=2)
                    plt.plot(cycle*0.1+Tstress_accum+(pulse-pulse_idx)*pulse_length[cycle] + np.arange(2*pulse_length[cycle]/(1.0+1.0+1.0), pulse_length[cycle]+0.0001, pulse_length[cycle]/(1.0+1.0+1.0)), 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)+1: data_points+1+(pulse-pulse_idx+1)*(1+1+1)], color=seg_color, linestyle='solid', linewidth = 0.4, alpha = 0.5)
                    

            Tstress_accum += write_time[cycle]*pulse_length[cycle]
            data_points += 1+write_time[cycle]*(1+1+1)
            pulse_idx += write_time[cycle]
            """

    """
            if(meas_dir == 'reversed'):
                print(last_IDSAT_above_threshold_round1[cycle])
                fig, ax = plt.subplots(nrows = 1, ncols = 1)
                n, bins, patches = ax.hist(last_IDSAT_above_threshold_round1[cycle], bins=40, normed = False)
                plt.ylabel('number of cells')
                plt.xlabel('last IDSAT before crossing threshold (uA)')
                plt.title('level-'+str(cycle+1)+', last IDSAT before crossing threshold (uA)', fontsize=8)
                plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_level'+str(cycle+1)+'_last_IDSAT_above_threshold.pdf')
                plt.close()

                print(delta_IDSAT_cross_threshold_round1[cycle])
                fig, ax = plt.subplots(nrows = 1, ncols = 1)
                n, bins, patches = ax.hist(delta_IDSAT_cross_threshold_round1[cycle], bins=40, normed = False)
                plt.ylabel('number of cells')
                plt.xlabel('delta IDSAT when crossing the threshold (uA)')
                plt.title('level-'+str(cycle+1)+', delta IDSAT during the pulse when crossing threshold (uA)', fontsize=8)
                plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_level'+str(cycle+1)+'_delta_IDSAT_cross_threshold.pdf')
                plt.close()
                """

        
    """
        #plt.xticks(t, t_label, rotation=30, fontsize=9)
        #plt.xticks(t, t_label, fontsize = 5)
        if (Imin > 0.01) and (Imax > 0.01):
            plt.ylim(Imin, Imax)
        if(t_range == []):
            plt.xlim(0, Tstress_accum+0.1*(PulseCycle-1))
        else:
            plt.xlim(t_range[0], t_range[1])
        plt.grid()
        #plt.legend([HCI_fig, PBTI_fig], ['HCI', 'PBTI: VDS=0'], loc = 'best')
        plt.xlabel('time (sec)')
        plt.ylabel('IDSAT (uA)')

        if(len(row_idx) == len(marker_lst)):
            plt.legend(row_axe, [str(row_num) for row_num in row_idx], loc='best', fontsize = 7, numpoints = 1)
                
        ##plt.subplots_adjust(bottom=0.15)
        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_'+meas_dir+title+'.pdf')
        #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_'+meas_dir+title+'_curve_fit_log.pdf')
        #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_'+meas_dir+title+'_curve_fit_log_bounded.pdf')
        #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_'+meas_dir+title+'_curve_fit_power.pdf')
        #plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_'+meas_dir+title+'_curve_fit_exp.pdf')
        ##plt.savefig(path_plot+'IDSAT_'+str(figN).zfill(3)+'_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+meas_dir+'.pdf')
        plt.close()
        """
        


    """
    if(Nrow == len(row_idx)):
        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+'\nrv1 2-bit MLC, Rows_remain & Rows_remain_round2 vs pulse, VGS='+str(VG)+'V, VDS='+str(VD)+'V', fontsize=10)
        pulse_idx = 0
        for cycle in np.arange(PulseCycle):
            round1, = plt.plot(pulse_idx+np.arange(0, write_time[cycle]+1), np.append(np.array([Nrow]), Rows_remain[pulse_idx:pulse_idx+write_time[cycle]]), marker='.', color = 'r')
            round2, = plt.plot(pulse_idx+np.arange(0, write_time[cycle]+1), np.append(np.array([Nrow]), Rows_remain_round2[pulse_idx:pulse_idx+write_time[cycle]]), marker='.', color = 'y')
            pulse_idx += write_time[cycle]

        t = np.arange(1, Total_PulseCycle+1)
        t_label = []
        for cycle in np.arange(PulseCycle):
            for pulse in np.arange(1, write_time[cycle]+1):
                t_label.append(str(pulse))

        plt.xticks(t, t_label, fontsize = 5)
        plt.ylim(0, Nrow)
        plt.xlim(0, Total_PulseCycle)
        plt.legend([round1, round2], ['round1', 'round2'])
        plt.grid()
        #plt.legend([HCI_fig, PBTI_fig], ['HCI', 'PBTI: VDS=0'], loc = 'best')
        plt.xlabel('pulse cycle')
        plt.ylabel('Number of rows remain')
        #plt.subplots_adjust(bottom=0.15)
        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+title+'_Rows_remain.pdf')
        #plt.savefig(path_plot+'IDSAT_'+str(figN).zfill(3)+'_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+meas_dir+'.pdf')
        plt.close()
        """

if __name__ == '__main__':
  main()

