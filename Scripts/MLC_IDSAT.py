
import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d

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
                plt.plot(cycle*0.4+np.append(np.arange(Tstress_accum-pulse_length[cycle], Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.0001, pulse_length[cycle]), np.array([Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.1])), 1e6*IDSAT[row][data_points: data_points+(1+1+write_time[cycle]+1)], color='r', linestyle='solid', marker='.')

            Tstress_accum += write_time[cycle]*pulse_length[cycle]
            data_points += 1+1+write_time[cycle]+1

        #plt.xticks(t, t_label, rotation=30, fontsize=9)
        plt.xticks(t, t_label)
        plt.ylim(Imin, Imax)
        plt.xlim(-pulse_length[0], t[-1]+0.1)
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
                plt.plot(cycle*0.4+np.arange(Tstress_accum+pulse_length[cycle], Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.0001, pulse_length[cycle]), 1e6*ID_prog[row][data_points: data_points+write_time[cycle]], color='r', linestyle='solid', marker='.')

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
                plt.plot(cycle*0.4+np.arange(Tstress_accum+pulse_length[cycle], Tstress_accum+write_time[cycle]*pulse_length[cycle]+0.0001, pulse_length[cycle]), 1e6*Isub[row][data_points: data_points+write_time[cycle]], color='g', linestyle='solid', marker='.')

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
        
        """
        #if (Imin < 0.01) and (Imax < 0.01):
            Imin = 1e6*np.amin(IDSAT[row_idx])
            Imax = 1e6*np.amax(IDSAT[row_idx])

        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+'\n, naive 2-bit MLC, IDSAT evolution, VGS='+str(VG)+'V, VDS='+str(VD)+'V\nIDSAT measured at VGS=0.8, VDS=0.8, '+meas_dir+' direction', fontsize=10)
        Tstress_accum = 0
        data_points = 0
        pulse_idx = 0
        for cycle in np.arange(PulseCycle):

            row_axe=[]
            marker_lst = ['.', 'o', 'D', '^', 'x', '+', 's', '*']

            for row in row_idx:
                if(len(row_idx) == len(marker_lst)):
                    row_marker = marker_lst[row - row_idx[0]]
                else:
                    row_marker = '.'

                axe, = plt.plot([cycle*0.1+Tstress_accum], 1e6*IDSAT[row][data_points], color='b', marker = row_marker)
                if(len(row_idx) == len(marker_lst)):
                    row_axe.append(axe)
                
                for pulse in np.arange(pulse_idx, pulse_idx+write_time[cycle]): 
                    if (Apply_Pulse[row][pulse] == 1):
                        seg_color = 'r'
                    if (Apply_Pulse[row][pulse] == 0):
                        seg_color = 'b'
                    plt.plot(cycle*0.1+Tstress_accum+(pulse-pulse_idx)*pulse_length[cycle] + np.array([0, pulse_length[cycle]/(1.0+1.0+1.0)]), 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)-1: data_points+1+(pulse-pulse_idx)*(1+1+1)+1], color=seg_color, linestyle='solid', linewidth = 0.4, alpha = 0.5)
                    plt.plot(cycle*0.1+Tstress_accum+(pulse-pulse_idx)*pulse_length[cycle] + np.arange(pulse_length[cycle]/(1.0+1.0+1.0), 2*pulse_length[cycle]/(1.0+1.0+1.0)+0.0001, pulse_length[cycle]/(1.0+1.0+1.0)), 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1): data_points+1+(pulse-pulse_idx+1)*(1+1+1)-1], color=seg_color, linestyle='solid', marker=row_marker, linewidth = 0.4, markersize=2)
                    plt.plot(cycle*0.1+Tstress_accum+(pulse-pulse_idx)*pulse_length[cycle] + np.arange(2*pulse_length[cycle]/(1.0+1.0+1.0), pulse_length[cycle]+0.0001, pulse_length[cycle]/(1.0+1.0+1.0)), 1e6*IDSAT[row][data_points+1+(pulse-pulse_idx)*(1+1+1)+1: data_points+1+(pulse-pulse_idx+1)*(1+1+1)], color=seg_color, linestyle='solid', linewidth = 0.4, alpha = 0.5)

            Tstress_accum += write_time[cycle]*pulse_length[cycle]
            data_points += 1+write_time[cycle]*(1+1+1)
            pulse_idx += write_time[cycle]

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
        
        #if (Imin < 0.01) and (Imax < 0.01):
            Imin = 1e6*np.amin(IDSAT[row_idx])
            Imax = 1e6*np.amax(IDSAT[row_idx])

        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+'\n, rv1 2-bit MLC, IDSAT evolution, VGS='+str(VG)+'V, VDS='+str(VD)+'V\nIDSAT measured at VGS=0.8, VDS=0.8, '+meas_dir+' direction', fontsize=10)
        Tstress_accum = 0
        data_points = 0
        pulse_idx = 0
        for cycle in np.arange(PulseCycle):

            row_axe=[]
            marker_lst = ['.', 'o', 'D', '^', 'x', '+', 's', '*']

            for row in row_idx:
                if(len(row_idx) == len(marker_lst)):
                    row_marker = marker_lst[row - row_idx[0]]
                else:
                    row_marker = '.'

                axe, = plt.plot([cycle*0.1+Tstress_accum], 1e6*IDSAT[row][data_points], color='b', marker = row_marker)
                if(len(row_idx) == len(marker_lst)):
                    row_axe.append(axe)
                
                for pulse in np.arange(pulse_idx, pulse_idx+write_time[cycle]): 
                    if ((round_idx[row][pulse] == 1) and (Next_Pulse[row][pulse] == 1)):
                        seg_color = 'r'
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

        #plt.xticks(t, t_label, rotation=30, fontsize=9)
        #plt.xticks(t, t_label, fontsize = 5)
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
                
        #plt.subplots_adjust(bottom=0.15)
        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+VGVD_char+'_'+meas_dir+title+'.pdf')
        #plt.savefig(path_plot+'IDSAT_'+str(figN).zfill(3)+'_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+meas_dir+'.pdf')
        plt.close()

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

