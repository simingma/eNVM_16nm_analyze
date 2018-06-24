
import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d

def RTN_ID_vs_time(chip, col, L, Nfin, VT_flavor, VG, VD, VGVD_char, Nrow, meas_direction_char, 
        ID_samples, sample_interval, sample_interval_char,
        path_data, path_plot, title='', Imin=0, Imax=0, tmin=0, tmax=0):

    """ Plot RTN in ID vs sample times.
    turning on WL=VG, setting VD, use DMM34410 ExtTrig = #ID_samples, with sample_interval, and integration time PLC
    tracking ID for a total time = #ID_samples x sample_interval.
    sample_frequency fs = 1/sample_interval (should be large enough as respect to RTN's power density spectrum: fs/2 >= PDS BandWidth
    """

    if os.path.isdir(path_plot) == False:
        os.mkdir(path_plot)

    ID = np.zeros((Nrow, ID_samples)) 
    f = open(path_data,'rU')
    for row in np.arange(0, Nrow):
        f.seek(0)
        ID_all = []
        ID_all.append(re.findall(r'WL\['+str(row)+']_ID=(-*\d*\.\d+)',f.read()))

        if (len(ID_all[0]) != ID_samples):
            print('data file error!\ngrabed '+str(len(ID_all[0]))+' WL['+str(row)+'] ID,\nbut expecting '+str(ID_samples)+' data\n')

        for sample in np.arange(0, ID_samples):
            ID[row][sample] = np.float64(ID_all[0][sample])
    
        if (Imin < 0.01) and (Imax < 0.01):
            Imin = 1e6*np.amin(ID)
            Imax = 1e6*np.amax(ID)
        if (tmin < 1e-7) and (tmax < 1e-7):
            tmin = 0
            tmax = ID_samples * sample_interval

        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+'\nID vs sample time, VGS='+str(VG)+'V, VDS='+str(VD)+'V, '+meas_direction_char+'\nsample interval = '+str(sample_interval)+'s, '+str(ID_samples)+' samples, chip'+str(chip)+', col['+str(col)+', row['+str(row)+']', fontsize=10)
        plt.plot(np.arange(0, ID_samples*sample_interval - sample_interval*0.01, sample_interval), 1e6*ID[row])
        #plt.xticks(t, t_label, rotation=30, fontsize=9)
        #plt.ylim(Imin, Imax)
        #plt.xlim(0, t[-1])
        plt.xlim(tmin, tmax)
        plt.grid()
        plt.xlabel('time (sec)')
        plt.ylabel('ID (uA)')
        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_row'+str(row)+'_ID_RTN_sample_interval_'+sample_interval_char+'_'+str(ID_samples)+'samples'+meas_direction_char+VGVD_char+title+'.pdf')
        plt.close()
    f.close()

if __name__ == '__main__':
  main()

