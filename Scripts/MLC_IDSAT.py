
import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d

def MLC_IDSAT(chip, col, L, Nfin, VT_flavor, VG, VD, Nrow, row_idx, 
        write_time, pulse_length, PulseCycle,  
        t, t_label,
        data_files, path_plot, VGVD_char, title, Imin=0, Imax=0):

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

if __name__ == '__main__':
  main()

