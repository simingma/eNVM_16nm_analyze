#!/usr/bin/python26 -tt

import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d

def Current_Component(chip, col, L, Nfin, VT_flavor, data_file, data_file_MUX_OFF, path_plot, VA_VB_direction, title, VB, VD_max, VD_step = 0.1, VG = 0, VS = 0):

    if os.path.isdir(path_plot) == False:
        os.mkdir(path_plot)

    VD_sweep = np.arange(0, VD_max + 0.0001, VD_step)

    if (data_file_MUX_OFF != ''):
        I_Drain_0 = []
        I_Psub_0 = []
        I_Source_0 = []
        f_0 = open(data_file_MUX_OFF, 'rU')
        I_Drain_0.append(re.findall(r'ID=(-*\d+\.\d+)',f_0.read()))
        f_0.seek(0, 0)
        I_Psub_0.append(re.findall(r'Isub=(-*\d+\.\d+)',f_0.read()))
        f_0.seek(0, 0)
        I_Source_0.append(re.findall(r'I_VDD_WL=(-*\d+\.\d+)',f_0.read()))
        #I_Source_0.append(re.findall(r'IS=(-*\d+\.\d+)',f_0.read()))
        f_0.close()
        print(len(I_Drain_0[0]), len(I_Psub_0[0]), len(I_Source_0[0])) #should all be len(VD) = (2.4/0.1+1)=25
    I_Drain = []
    I_Psub = []
    I_Source = []
    f = open(data_file, 'rU')
    I_Drain.append(re.findall(r'ID=(-*\d+\.\d+)',f.read()))
    f.seek(0, 0)
    I_Psub.append(re.findall(r'Isub=(-*\d+\.\d+)',f.read()))
    f.seek(0, 0)
    I_Source.append(re.findall(r'I_VDD_WL=(-*\d+\.\d+)',f.read()))
    #I_Source.append(re.findall(r'IS=(-*\d+\.\d+)',f.read()))
    f.close()
    print(len(I_Drain[0]), len(I_Psub[0]), len(I_Source[0])) #should all be len(VD) = (2.4/0.1+1)=25

    Isense_Drain = np.zeros(len(VD_sweep))
    Isense_Psub = np.zeros(len(VD_sweep))
    Isense_Source = np.zeros(len(VD_sweep))
    if (data_file_MUX_OFF != ''):
        Isense_Drain_0 = np.zeros(len(VD_sweep))
        Isense_Psub_0 = np.zeros(len(VD_sweep))
        Isense_Source_0 = np.zeros(len(VD_sweep))
    for k in np.arange(0, len(VD_sweep)):
        Isense_Drain[k] = np.float64(I_Drain[0][k])
        Isense_Psub[k] = np.float64(I_Psub[0][k])
        Isense_Source[k] = np.float64(I_Source[0][k])
        if (data_file_MUX_OFF != ''):
            Isense_Drain_0[k] = np.float64(I_Drain_0[0][k])
            Isense_Psub_0[k] = np.float64(I_Psub_0[0][k])
            Isense_Source_0[k] = np.float64(I_Source_0[0][k])

    ID, = plt.plot(VD_sweep, 1e6*Isense_Drain, color = 'r', linestyle = 'dashed', marker = '.')
    Isub, = plt.plot(VD_sweep, 1e6*Isense_Psub, color = 'g', linestyle = 'dashed', marker = '.')
    IS, = plt.plot(VD_sweep, 1e6*Isense_Source, color = 'y', linestyle = 'dashed', marker = '.')
    Itot, = plt.plot(VD_sweep, 1e6*(Isense_Drain + Isense_Source - Isense_Psub), color = 'm')
    if (data_file_MUX_OFF != ''):
        ID, = plt.plot(VD_sweep, 1e6*Isense_Drain_0, color = 'r', linestyle = 'dotted', marker = '+')
        Isub, = plt.plot(VD_sweep, 1e6*Isense_Psub_0, color = 'g', linestyle = 'dotted', marker = '+')
        IS, = plt.plot(VD_sweep, 1e6*Isense_Source_0, color = 'y', linestyle = 'dotted', marker = '+')
        ID, = plt.plot(VD_sweep, 1e6*(Isense_Drain-Isense_Drain_0), color = 'r')
        Isub, = plt.plot(VD_sweep, 1e6*(Isense_Psub-Isense_Psub_0), color = 'g')
        IS, = plt.plot(VD_sweep, 1e6*(Isense_Source-Isense_Source_0), color = 'y')
    plt.legend([ID, Isub, IS, Itot], ['ID', 'Isub', 'I_VDD_WL', 'ID+Isub+I_VDD_WL'], loc = 'best')
    #plt.legend([ID, Isub, IS], ['ID', 'Isub', 'IS'], loc = 'best')
    plt.xlabel('VD (V)')
    plt.ylabel('Current (uA)')
    if (data_file_MUX_OFF != ''):
        plt.title('BJT (S-B-D) current, VS='+str(VS)+', VG='+str(VG)+', VB='+str(VB)+', sweeping VD from 0 to '+str(VD_max)+'\nIS IB(i.e. Isub) ID with background (MUX-OFF) subtracted', fontsize = 10)
    else:
        plt.title('BJT (S-B-D) current, VS='+str(VS)+', VG='+str(VG)+', VB='+str(VB)+', sweeping VD from 0 to '+str(VD_max)+', IS IB(i.e. Isub) ID', fontsize = 10)
    #plt.title(fname, fontsize = 12)
    plt.grid()
    #plt.xlim(xmin = VD_min)
    #plt.ylim(ymin = -5)
    #plt.savefig(path_plot+str(figN).zfill(3)+fname+'_LVT.pdf')
    plt.savefig(path_plot+'BJT-current_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+title+VA_VB_direction+'.pdf')
    #plt.semilogy(VAB, 1e6*I_leak_Psub_VAsource_VBdrain , color='green', linestyle='solid')
    plt.close()

def main():
    chip = 17
    col = 33
    L = 16
    Nfin = 2
    VT_flavor = 'ULVT'
    path_plot = '../Plots/chip17_BJT_current/'
    #Current_Component(chip, col, L, Nfin, VT_flavor, '../Data/chip17/Current_Components_VDD_IO_2p4_Vg0_VsVb0_Chip17_Col33_VAsource_VBdrain', '../Data/chip17/Current_Components_VDD_IO_2p4_Vg0_Vs0_Vb0_Chip17_Col33_MUX-OFF', path_plot, '_VAsource_VBdrain', '_Vg0_Vs0_Vb0_sweepVd2p4_VDD_WL_0p8_VDD_IO_2p4', 0, 2.4) 
    #Current_Component(chip, col, L, Nfin, VT_flavor, '../Data/chip17/Current_Components_VDD_IO_2p4_Vg0_Vs0_Vb1p4_Chip17_Col33_VAsource_VBdrain_1min-VD2p4', '../Data/chip17/Current_Components_VDD_IO_2p4_Vg0_Vs0_Vb1p4_Chip17_Col33_MUX-OFF', path_plot, '_VAsource_VBdrain', '_Vg0_Vs0_Vb1p4_sweepVd2p4_VDD_WL_0p8_VDD_IO_2p4', 1.4, 2.4) 
    #Current_Component(chip, col, L, Nfin, VT_flavor, '../Data/chip17/Current_Components_VDD_IO_2p4_Vg0_Vs0_Vb1p4_VDD_WL-UnPlug_Chip17_Col33_VAsource_VBdrain_SweepVD2p4', '../Data/chip17/Current_Components_VDD_IO_2p4_Vg0_Vs0_Vb1p4_VDD_WL-UnPlug_Chip17_Col33_MUX-OFF', path_plot, '_VAsource_VBdrain', '_Vg0_Vs0_Vb1p4_sweepVd2p4_VDD_WL-UnPlug_VDD_IO_2p4', 1.4, 2.4) 
    #Current_Component(chip, col, L, Nfin, VT_flavor, '../Data/chip17/Current_Components_VDD_IO_2p4_Vg0_Vs0_Vb1p4_VDD_WL-1p4_Chip17_Col33_VAsource_VBdrain_SweepVD2p4', '', path_plot, '_VAsource_VBdrain', '_Vg0_Vs0_Vb1p4_sweepVd2p4_VDD_WL-1p4_VDD_IO_2p4', 1.4, 2.4) 
    #Current_Component(chip, col, L, Nfin, VT_flavor, '../Data/chip17/Current_Components_VDD_IO_2p4_Vg0_Vs0_Vb1p4_VDD_WL-0p3_Chip17_Col33_VAsource_VBdrain_SweepVD2p4', '', path_plot, '_VAsource_VBdrain', '_Vg0_Vs0_Vb1p4_sweepVd2p4_VDD_WL-0p3_VDD_IO_2p4', 1.4, 2.4) 
    Current_Component(chip, col, L, Nfin, VT_flavor, '../Data/chip17/Current_Components_VDD_IO_2p4_Vg0_Vs0_Vb1p4_VDD_WL-0p8_Chip17_Col33_VAsource_VBdrain_I_VDD_WL', '', path_plot, '_VAsource_VBdrain', '_Vg0_Vs0_Vb1p4_sweepVd2p4_VDD_WL-0p8_VDD_IO_2p4_Measure-I_VDD_WL', 1.4, 2.4) 

if __name__ == '__main__':
  main()    
