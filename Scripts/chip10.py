
import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d
from VG_ConstPulse_IDSAT import IDSAT
from I_V_curves import IDS_VGS_stress, IDS_VGS, IDSAT_vs_row
from VG_ConstPulse_horizontal_hist import IDSAT_horizontal_hist
from VG_ConstPulse_separation import IDSAT_separation
from Charge_Pumping import Charge_Pumping, Charge_Pumping_compare
from MLC_IDSAT import MLC_IDSAT_characterization, MLC_IDSAT_algorithm_naivete

def main():
    
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

    MLC_IDSAT_algorithm_naivete(10, 21, 36, 2, 'ULVT', 1.8, 2.4, 128, range(0, 128), [16, 18, 15], [0.01, 0.04, 0.2], 3, [0, 0.16, 0.26, 0.98, 1.08, 4.08], ['0', '0.16', '0.16', '0.88', '0.88', '3.88'], ['../Data/chip10/MLC_programming_Chip10_Col21_10msPULSE_VG1p8_VD2p4_VAsource_VBdrain_01', '../Data/chip10/MLC_programming_Chip10_Col21_40msPULSE_VG1p8_VD2p4_VAsource_VBdrain_02', '../Data/chip10/MLC_programming_Chip10_Col21_200msPULSE_VG1p8_VD2p4_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD2p4', '_Naivete_cycle010203')

    #t_label = []
    #for t in np.arange(0, 0.01*16 + 0.0001, 0.01):
    #    t_label.append(str(t))
    ##MLC_IDSAT_algorithm_naivete(10, 21, 36, 2, 'ULVT', 1.8, 2.4, 128, range(0, 128), [16], [0.01], 1, np.arange(0, 0.01*16+0.0001, 0.01), t_label, ['../Data/chip10/MLC_programming_Chip10_Col21_10msPULSE_VG1p8_VD2p4_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD2p4', '_Naivete_cycle01')
    #for row_start in np.arange(0, 128, 8):
    #    MLC_IDSAT_algorithm_naivete(10, 21, 36, 2, 'ULVT', 1.8, 2.4, 128, range(row_start, row_start+8), [16], [0.01], 1, np.arange(0, 0.01*16+0.0001, 0.01), t_label, ['../Data/chip10/MLC_programming_Chip10_Col21_10msPULSE_VG1p8_VD2p4_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD2p4', '_Naivete_cycle01_row'+str(row_start)+'_to_'+str(row_start+7))


    #t_label = []
    #for t in np.arange(0, 0.04*18 + 0.0001, 0.04):
    #    t_label.append(str(0.01*16 + t))
    ##MLC_IDSAT_algorithm_naivete(10, 21, 36, 2, 'ULVT', 1.8, 2.4, 128, range(0, 128), [18], [0.04], 1, np.arange(0, 0.04*18+0.0001, 0.04), t_label, ['../Data/chip10/MLC_programming_Chip10_Col21_40msPULSE_VG1p8_VD2p4_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD2p4', '_Naivete_cycle02')
    #for row_start in np.arange(0, 128, 8):
    #    MLC_IDSAT_algorithm_naivete(10, 21, 36, 2, 'ULVT', 1.8, 2.4, 128, range(row_start, row_start+8), [18], [0.04], 1, np.arange(0, 0.04*18+0.0001, 0.04), t_label, ['../Data/chip10/MLC_programming_Chip10_Col21_40msPULSE_VG1p8_VD2p4_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD2p4', '_Naivete_cycle02_row'+str(row_start)+'_to_'+str(row_start+7))


    #t_label = []
    #for t in np.arange(0, 0.2*15 + 0.0001, 0.2):
    #    t_label.append(str(0.01*16 + 0.04*18 + t))
    ##MLC_IDSAT_algorithm_naivete(10, 21, 36, 2, 'ULVT', 1.8, 2.4, 128, range(0, 128), [15], [0.2], 1, np.arange(0, 0.2*15+0.0001, 0.2), t_label, ['../Data/chip10/MLC_programming_Chip10_Col21_200msPULSE_VG1p8_VD2p4_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD2p4', '_Naivete_cycle03')
    #for row_start in np.arange(0, 128, 8):
    #    MLC_IDSAT_algorithm_naivete(10, 21, 36, 2, 'ULVT', 1.8, 2.4, 128, range(row_start, row_start+8), [15], [0.2], 1, np.arange(0, 0.2*15+0.0001, 0.2), t_label, ['../Data/chip10/MLC_programming_Chip10_Col21_200msPULSE_VG1p8_VD2p4_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD2p4', '_Naivete_cycle03_row'+str(row_start)+'_to_'+str(row_start+7))

    #MLC_IDSAT_characterization(10, 18, 36, 2, 'ULVT', 1.8, 2.0, 32, range(0, 16) , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip10/Chip10_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col18_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col18_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD2p0', '_cycle010203', 38, 112, 1)
    #MLC_IDSAT_characterization(10, 18, 36, 2, 'ULVT', 1.8, 2.4, 32, range(16, 32), [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip10/Chip10_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col18_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col18_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD2p4', '_cycle010203', 16, 110, 1)

    #MLC_IDSAT_characterization(10, 24, 20, 2, 'ULVT', 1.8, 1.8, 32, range(0, 16) , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip10/Chip10_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col24_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col24_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD1p8', '_cycle010203', 44, 133, 1)
    #MLC_IDSAT_characterization(10, 24, 20, 2, 'ULVT', 1.8, 2.2, 32, range(16, 32), [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip10/Chip10_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col24_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col24_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD2p2', '_cycle010203', 14, 133, 1)

    #MLC_IDSAT_characterization(10, 30, 16, 2, 'ULVT', 1.8, 1.7, 32, range(0, 16) , [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip10/Chip10_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col30_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col30_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD1p7', '_cycle010203', 50, 135, 1)
    #MLC_IDSAT_characterization(10, 30, 16, 2, 'ULVT', 1.8, 2.0, 32, range(16, 32), [40, 20, 12], [0.01, 0.04, 0.2], 3, [0, 0.4, 0.8, 1.6, 2.0, 2.6, 3.2, 3.8, 4.4], ['0', '0.4', '0.4', '1.2', '1.2', '1.8', '2.4', '3.0', '3.6'], ['../Data/chip10/Chip10_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col30_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02', '../Data/chip10/Chip10_Col30_HCI_12x200ms_stress_VG_ConstPulse_VAsource_VBdrain_03'], '../Plots/chip10/', 'VG1p8_VD2p0', '_cycle010203', 20, 140, 1)



    #MLC_IDSAT_characterization(10, 18, 36, 2, 'ULVT', 1.8, 2.0, 32, range(0, 16) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip10/Chip10_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD2p0', '_cycle01', 38, 112, 1)
    #MLC_IDSAT_characterization(10, 18, 36, 2, 'ULVT', 1.8, 2.4, 32, range(16, 32), [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip10/Chip10_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD2p4', '_cycle01', 16, 110, 1)

    #MLC_IDSAT_characterization(10, 24, 20, 2, 'ULVT', 1.8, 1.8, 32, range(0, 16) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip10/Chip10_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD1p8', '_cycle01', 44, 133, 1)
    #MLC_IDSAT_characterization(10, 24, 20, 2, 'ULVT', 1.8, 2.2, 32, range(16, 32), [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip10/Chip10_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD2p2', '_cycle01', 14, 133, 1)

    #MLC_IDSAT_characterization(10, 30, 16, 2, 'ULVT', 1.8, 1.7, 32, range(0, 16) , [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip10/Chip10_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD1p7', '_cycle01', 50, 135, 1)
    #MLC_IDSAT_characterization(10, 30, 16, 2, 'ULVT', 1.8, 2.0, 32, range(16, 32), [40], [0.01], 1, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5], ['0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', 'recover'], ['../Data/chip10/Chip10_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01'], '../Plots/chip10/', 'VG1p8_VD2p0', '_cycle01', 20, 140, 1)



    #MLC_IDSAT_characterization(10, 18, 36, 2, 'ULVT', 1.8, 2.0, 32, range(0, 16) , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip10/Chip10_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col18_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD2p0', '_cycle0102', 38, 112, 1)
    #MLC_IDSAT_characterization(10, 18, 36, 2, 'ULVT', 1.8, 2.4, 32, range(16, 32), [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip10/Chip10_Col18_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col18_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD2p4', '_cycle0102', 16, 110, 1)
                                                                                                            
    #MLC_IDSAT_characterization(10, 24, 20, 2, 'ULVT', 1.8, 1.8, 32, range(0, 16) , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip10/Chip10_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col24_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD1p8', '_cycle0102', 44, 133, 1)
    #MLC_IDSAT_characterization(10, 24, 20, 2, 'ULVT', 1.8, 2.2, 32, range(16, 32), [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip10/Chip10_Col24_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col24_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD2p2', '_cycle0102', 14, 133, 1)
                                                                                                            
    #MLC_IDSAT_characterization(10, 30, 16, 2, 'ULVT', 1.8, 1.7, 32, range(0, 16) , [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip10/Chip10_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col30_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD1p7', '_cycle0102', 50, 135, 1)
    #MLC_IDSAT_characterization(10, 30, 16, 2, 'ULVT', 1.8, 2.0, 32, range(16, 32), [40, 20], [0.01, 0.04], 2, [0, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.7], ['0', '0.2', '0.4', '0.4', '0.6', '0.8', '1.0', '1.2', 'recover'], ['../Data/chip10/Chip10_Col30_HCI_40x10ms_stress_VG_ConstPulse_VAsource_VBdrain_01', '../Data/chip10/Chip10_Col30_HCI_20x40ms_stress_VG_ConstPulse_VAsource_VBdrain_02'], '../Plots/chip10/', 'VG1p8_VD2p0', '_cycle0102', 20, 140, 1)


#    for col in range(36):
#        IDS_VGS(10, col, col_list[col][0], col_list[col][1], col_list[col][2], col_list[col][3], ['../Data/chip10/Fresh_Chip10_Col'+str(col).zfill(2)+'_Ids_Vgs_VAsource_VBdrain'], ['b'], '../Plots/chip10/', 'Fresh_Ids-Vgs_VaS-VbD_', range(0, col_list[col][3]), 'Fresh IDS-VGS, forward' , col_list[col][4])
#        IDS_VGS(10, col, col_list[col][0], col_list[col][1], col_list[col][2], col_list[col][3], ['../Data/chip10/Fresh_Chip10_Col'+str(col).zfill(2)+'_Ids_Vgs_VAdrain_VBsource'], ['b'], '../Plots/chip10/', 'Fresh_Ids-Vgs_VaD-VbS_', range(0, col_list[col][3]), 'Fresh IDS-VGS, reversed', col_list[col][4])

    #IDS_VGS(10, 21, 36, 2, 'ULVT', 128, ['../Data/chip10/Fresh_Chip10_Col21_Ids_Vgs_VAsource_VBdrain', '../Data/chip10/MLC_Chip10_Col21_10msPULSE_VG1p8_VD2p4_Ids_Vgs_VAsource_VBdrain_01'], ['b', 'y'], '../Plots/chip10/', 'Fresh_vs_MLC01_VG1p8_VD2p4_IDS-VGS_VaS-VbD_', range(0, 128), 'Fresh vs MLC-01 (VG=1.8, VD=2.4)\nMLC-01: 10ms WL pulse, IDSAT threshold = 80uA, forward' , col_list[21][4], ['fresh', 'MLC-01']) 
    #IDS_VGS(10, 21, 36, 2, 'ULVT', 128, ['../Data/chip10/Fresh_Chip10_Col21_Ids_Vgs_VAdrain_VBsource', '../Data/chip10/MLC_Chip10_Col21_10msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_01'], ['b', 'y'], '../Plots/chip10/', 'Fresh_vs_MLC01_VG1p8_VD2p4_IDS-VGS_VaD-VbS_', range(0, 128), 'Fresh vs MLC-01 (VG=1.8, VD=2.4)\nMLC-01: 10ms WL pulse, IDSAT threshold = 80uA, reversed', col_list[21][4], ['fresh', 'MLC-01']) 

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
    
