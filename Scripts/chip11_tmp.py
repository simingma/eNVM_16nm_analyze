
import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d
from VG_ConstPulse_IDSAT import IDSAT
from I_V_curves import hist_IDS_VGS, IDS_VGS_stress, IDS_VGS, IDSAT_vs_row
from VG_ConstPulse_horizontal_hist import IDSAT_horizontal_hist
from VG_ConstPulse_separation import IDSAT_separation
from Charge_Pumping import Charge_Pumping, Charge_Pumping_compare
from MLC_IDSAT import MLC_IDSAT_characterization, MLC_IDSAT_algorithm_naivete, MLC_IDSAT_algorithm_rv1

def main():
    hist_IDS_VGS(0, 11, 21, 36, 2, 'ULVT', 128, ['../Data/chip11/Fresh_Chip11_Col21_Ids_Vgs_VAdrain_VBsource', '../Data/chip11/MLC_Chip11_Col21_2msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_01', '../Data/chip11/MLC_Chip11_Col21_10msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_02', '../Data/chip11/MLC_Chip11_Col21_40msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_03', '../Data/chip11/MLC_Chip11_Col21_200msPULSE_VG1p8_VD2p4_Ids_Vgs_VAdrain_VBsource_04'], ['b', 'y', 'r', 'k', 'g'], '../Plots/chip11/', 'Hist-IDSAT_MLC-rv1-01020304_reverse-read_', range(0, 128), 'MLC programming (VGS=1.8, VDS=2.4), histogram of read-IDSAT (VGS=VDS=0.8V)', 0, 136, 0, 136, 1000)
    

if __name__ == '__main__':
  main()

