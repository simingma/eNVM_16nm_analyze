
import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d
from VG_ConstPulse_IDSAT import IDSAT
from I_V_curves import IDS_VGS_stress, IDS_VGS
from VG_ConstPulse_horizontal_hist import IDSAT_horizontal_hist
from VG_ConstPulse_separation import IDSAT_separation
from Charge_Pumping import Charge_Pumping, Charge_Pumping_compare

def main():

    IDS_VGS(9, 21, 36, 2, 'ULVT', 128, ['../Data/chip09/Fresh_Chip09_Col21_Ids_Vgs_VAsource_VBdrain', '../Data/chip09/HCI_VG1p8_VD2p0_1x1p2s_Chip09_Col21_Ids_Vgs_VAsource_VBdrain_01'], ['b', 'r'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_IDS-VGS_VaS-VbD_', range(0, 128), 'Fresh (b) vs Stress01 (VG=1.8, VD=2.0, 1.2s) (r), forward', 130) 
    IDS_VGS(9, 21, 36, 2, 'ULVT', 128, ['../Data/chip09/Fresh_Chip09_Col21_Ids_Vgs_VAdrain_VBsource', '../Data/chip09/HCI_VG1p8_VD2p0_1x1p2s_Chip09_Col21_Ids_Vgs_VAdrain_VBsource_01'], ['b', 'r'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_IDS-VGS_VaD-VbS_', range(0, 128), 'Fresh (b) vs Stress01 (VG=1.8, VD=2.0, 1.2s) (r), reversed', 130) 

    IDS_VGS(9, 23, 36, 2, 'SVT', 128, ['../Data/chip09/Fresh_Chip09_Col23_Ids_Vgs_VAsource_VBdrain', '../Data/chip09/HCI_VG1p8_VD2p4_1x1p2s_Chip09_Col23_Ids_Vgs_VAsource_VBdrain_01'], ['b', 'r'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_IDS-VGS_VaS-VbD_', range(0, 128), 'Fresh (b) vs Stress01 (VG=1.8, VD=2.4, 1.2s) (r), forward', 95) 
    IDS_VGS(9, 23, 36, 2, 'SVT', 128, ['../Data/chip09/Fresh_Chip09_Col23_Ids_Vgs_VAdrain_VBsource', '../Data/chip09/HCI_VG1p8_VD2p4_1x1p2s_Chip09_Col23_Ids_Vgs_VAdrain_VBsource_01'], ['b', 'r'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_IDS-VGS_VaD-VbS_', range(0, 128), 'Fresh (b) vs Stress01 (VG=1.8, VD=2.4, 1.2s) (r), reversed', 95) 

    IDS_VGS(9, 33, 16, 2, 'ULVT', 128, ['../Data/chip09/Fresh_Chip09_Col33_Ids_Vgs_VAsource_VBdrain', '../Data/chip09/HCI_VG1p8_VD2p0_1x1p2s_Chip09_Col33_Ids_Vgs_VAsource_VBdrain_01'], ['b', 'r'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_IDS-VGS_VaS-VbD_', range(0, 128), 'Fresh (b) vs Stress01 (VG=1.8, VD=2.0, 1.2s) (r), forward', 170) 
    IDS_VGS(9, 33, 16, 2, 'ULVT', 128, ['../Data/chip09/Fresh_Chip09_Col33_Ids_Vgs_VAdrain_VBsource', '../Data/chip09/HCI_VG1p8_VD2p0_1x1p2s_Chip09_Col33_Ids_Vgs_VAdrain_VBsource_01'], ['b', 'r'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_IDS-VGS_VaD-VbS_', range(0, 128), 'Fresh (b) vs Stress01 (VG=1.8, VD=2.0, 1.2s) (r), reversed', 170) 

    IDS_VGS(9, 35, 16, 2, 'SVT', 128, ['../Data/chip09/Fresh_Chip09_Col35_Ids_Vgs_VAsource_VBdrain', '../Data/chip09/HCI_VG1p8_VD1p7_1x1p2s_Chip09_Col35_Ids_Vgs_VAsource_VBdrain_01'], ['b', 'r'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_IDS-VGS_VaS-VbD_', range(0, 128), 'Fresh (b) vs Stress01 (VG=1.8, VD=1.7, 1.2s) (r), forward', 125) 
    IDS_VGS(9, 35, 16, 2, 'SVT', 128, ['../Data/chip09/Fresh_Chip09_Col35_Ids_Vgs_VAdrain_VBsource', '../Data/chip09/HCI_VG1p8_VD1p7_1x1p2s_Chip09_Col35_Ids_Vgs_VAdrain_VBsource_01'], ['b', 'r'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_IDS-VGS_VaD-VbS_', range(0, 128), 'Fresh (b) vs Stress01 (VG=1.8, VD=1.7, 1.2s) (r), reversed', 125) 

    #Charge_Pumping_compare([9, 9], ['Fresh', 'VG=1.8, VD=2.0, 1x1.2sec'], 21, 36, 2, 'ULVT', 128, ['../Data/chip09/Fresh_Chip09_Col21_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM', '../Data/chip09/HCIstress01_1x1p2s_Chip09_Col21_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_VG1p8_VD2p0_1p2s_Charge_Pumping_VSS_WL_0_VDD_WL_1p6', 'Fresh vs Stress (VG=1.8, VD=2.0, 1x1.2s), col[21], 5MHz Charge Pumping\n', -1.6, 0, 0.1, 17, 0, 1.6, 60, [2, 2], [[5000000, 1000], [5000000, 1000]], 0, [1, 1], '1kHz', [[0], [0]])

    #Charge_Pumping_compare([9, 9], ['Fresh', 'VG=1.8, VD=2.4, 1x1.2sec'], 23, 36, 2, 'SVT', 128, ['../Data/chip09/Fresh_Chip09_Col23_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM', '../Data/chip09/HCIstress01_1x1p2s_Chip09_Col23_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_VG1p8_VD2p4_1p2s_Charge_Pumping_VSS_WL_0_VDD_WL_1p6', 'Fresh vs Stress (VG=1.8, VD=2.4, 1x1.2s), col[23], 5MHz Charge Pumping\n', -1.6, 0, 0.1, 17, 0, 1.6, 60, [2, 2], [[5000000, 1000], [5000000, 1000]], 0, [1, 1], '1kHz', [[0], [0]])

    #Charge_Pumping_compare([9, 9], ['Fresh', 'VG=1.8, VD=2.0, 1x1.2sec'], 33, 16, 2, 'ULVT', 128, ['../Data/chip09/Fresh_Chip09_Col33_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM', '../Data/chip09/HCIstress01_1x1p2s_Chip09_Col33_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_VG1p8_VD2p0_1p2s_Charge_Pumping_VSS_WL_0_VDD_WL_1p6', 'Fresh vs Stress (VG=1.8, VD=2.0, 1x1.2s), col[33], 5MHz Charge Pumping\n', -1.6, 0, 0.1, 17, 0, 1.6, 60, [2, 2], [[5000000, 1000], [5000000, 1000]], 0, [1, 1], '1kHz', [[0], [0]])

    #Charge_Pumping_compare([9, 9], ['Fresh', 'VG=1.8, VD=1.7, 1x1.2sec'], 35, 16, 2, 'SVT', 128, ['../Data/chip09/Fresh_Chip09_Col35_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM', '../Data/chip09/HCIstress01_1x1p2s_Chip09_Col35_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM'], '../Plots/chip09/', 'Fresh_vs_HCIstress01_VG1p8_VD1p7_1p2s_Charge_Pumping_VSS_WL_0_VDD_WL_1p6', 'Fresh vs Stress (VG=1.8, VD=1.7, 1x1.2s), col[35], 5MHz Charge Pumping\n', -1.6, 0, 0.1, 17, 0, 1.6, 60, [2, 2], [[5000000, 1000], [5000000, 1000]], 0, [1, 1], '1kHz', [[0], [0]])

    #Charge_Pumping(8, 21, 36, 2, 'ULVT', 128, '../Data/chip09/Fresh_Chip09_Col21_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM', '../Plots/chip09/', 'Fresh_Charge_Pumping_VSS_WL_0_VDD_WL_1p6', '', -1.6, 0, 0.1, 17, 0, 1.6, 60, 2, [5000000, 1000], ['5MHz', '1kHz'], 0, [0, 1], 1, [0], 0)

    #Charge_Pumping(8, 23, 36, 2, 'SVT', 128, '../Data/chip09/Fresh_Chip09_Col23_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM', '../Plots/chip09/', 'Fresh_Charge_Pumping_VSS_WL_0_VDD_WL_1p6', '', -1.6, 0, 0.1, 17, 0, 1.6, 60, 2, [5000000, 1000], ['5MHz', '1kHz'], 0, [0, 1], 1, [0], 0)

    #Charge_Pumping(8, 33, 16, 2, 'ULVT', 128, '../Data/chip09/Fresh_Chip09_Col33_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM', '../Plots/chip09/', 'Fresh_Charge_Pumping_VSS_WL_0_VDD_WL_1p6', '', -1.6, 0, 0.1, 17, 0, 1.6, 60, 2, [5000000, 1000], ['5MHz', '1kHz'], 0, [0, 1], 1, [0], 0)

    #Charge_Pumping(8, 35, 16, 2, 'SVT', 128, '../Data/chip09/Fresh_Chip09_Col35_60Pumping_SweepVSVBVD_VSS_WL_0_VDD_WL_1p6_ELTM', '../Plots/chip09/', 'Fresh_Charge_Pumping_VSS_WL_0_VDD_WL_1p6', '', -1.6, 0, 0.1, 17, 0, 1.6, 60, 2, [5000000, 1000], ['5MHz', '1kHz'], 0, [0, 1], 1, [0], 0)

if __name__ == '__main__':
  main()
    
