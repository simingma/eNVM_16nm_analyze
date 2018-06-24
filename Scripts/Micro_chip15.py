
import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from Micro_I_V_curves import hist_IDS_VGS, IDS_VGS

def main():
    
    IDS_VGS(15, 33, 16, 2, 'ULVT', 128, ['/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/Fresh_Chip15_Col33_Ids_Vgs_VAdrain_VBsource', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_1msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_01', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_3msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_02', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_9msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_03', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_27msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_04', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_81msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_05', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_243msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_06', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_729msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_07'], ['b', 'y', 'r', 'k', 'g', 'm', 'navy', 'blueviolet'], 'your_save_directory', 'test_Fresh_vs_MLC-1-7_VG1p8_VD2p0_IDS-VGS_VaD-VbS_', range(0, 128), 'Fresh vs MLC-1-7 (VG=1.8, VD=2.0)\nMLC-{1, 3, 9, 27, 81, 243, 729}ms WL pulses, IDSAT threshold = {95, 80, 65, 50, 35, 20, 5}uA, reversed', 165, []) 


    hist_IDS_VGS(0, 15, 33, 16, 2, 'ULVT', 128, ['/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/Fresh_Chip15_Col33_Ids_Vgs_VAdrain_VBsource', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_1msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_01', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_3msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_02', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_9msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_03', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_27msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_04', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_81msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_05', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_243msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_06', '/group/vlsiarch/siming/eNVM_16nm_analyze/Data/chip15/MLC_Chip15_Col33_729msPULSE_VG1p8_VD2p0_Ids_Vgs_VAdrain_VBsource_07'], ['b', 'y', 'r', 'k', 'g', 'm', 'navy', 'blueviolet'], 'your_save_directory', 'test_Hist-IDSAT_MLC-1-7_reverse-read_', range(0, 128), 'MLC programming {1, 3, 9, 27, 81, 243, 729}ms pulses, VGS=1.8, VDS=2.0 for level[1-7]\nhistogram of read-IDSAT (VGS=VDS=0.8V)', 0, 165, 0, 165, 1000)


if __name__ == '__main__':
  main()
    
