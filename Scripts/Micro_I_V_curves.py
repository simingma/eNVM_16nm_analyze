
import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d


def hist_IDS_VGS(VG_idx, chip, col, L, Nfin, VT_flavor, Nrow, data_files, colors, path_plot, plot_file, row_idx, title, Ymin, Ymax, Imin, Imax, Num_bins = 50, VD = 0.8, VG_min = 0.2, VG_max = 0.8):

    if os.path.isdir(path_plot) == False:
        os.mkdir(path_plot)

    VDD_WL = np.arange(VG_max, VG_min -0.0001, -0.05)
    #Ymax = 0
#    for direction in ['VAsource_VBdrain', 'VAdrain_VBsource']:
    fig_h, ax_h = plt.subplots(nrows = 1, ncols = 1)
    n_max = 0
    #mean = np.zeros(len(data_files))
    #std_dev = np.zeros(len(data_files))
    for (data_file, color) in zip(data_files, colors):
        Isense = []
        #f=open(path_data+'Fresh_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_Ids_Vgs_'+direction,'rU')
        f=open(data_file,'rU')
        Isense.append(re.findall(r'Isense=(-*\d*\.\d+)',f.read()))
        f.close()
        print(len(Isense), len(Isense[0]))
        if (len(Isense) != 1 or len(Isense[0]) != (Nrow+1)*len(VDD_WL)):
            print('data file error!')
        #I_leak_VAsource_VBdrain = np.zeros(len(VDD_WL))
        #I_leak_VAdrain_VBsource = np.zeros(len(VDD_WL))
        #for k in np.arange(len(VDD_WL)):
        #    I_leak_VAsource_VBdrain[k] = np.float64(Isense_VAsource_VBdrain[k][0])
        #    I_leak_VAdrain_VBsource[k] = np.float64(Isense_VAdrain_VBsource[k][0])
        IDS_VGS = np.zeros((Nrow, len(VDD_WL)))
        for row in np.arange(0, Nrow):
            IDS_VGS[row][:] = np.array(np.float64(Isense[0][(row+1)*len(VDD_WL): (row+2)*len(VDD_WL)]))

        #Ymax = max(1e6 * np.amax(IDS_VGS), Ymax)
        #plt.figure(figN)
        #plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', #'+str(Nrow)+' fresh nmos IDS-VGS curves', fontsize=10)

        #if Imin==0 or Imax==0:
        #    Imin=np.amin(1e6*IDS_VGS[row_idx, [VG_idx]*len(row_idx)])
        #    Imax=np.amax(1e6*IDS_VGS[row_idx, [VG_idx]*len(row_idx)])
        n, bins, patches = ax_h.hist(1e6*IDS_VGS[row_idx, [VG_idx]*len(row_idx)], bins=Num_bins, normed = False, range=(Imin, Imax), orientation = 'horizontal', color = color, edgecolor='none')
        mean = np.mean(1e6*IDS_VGS[row_idx, [VG_idx]*len(row_idx)])
        std_dev = np.std(1e6*IDS_VGS[row_idx, [VG_idx]*len(row_idx)])
        print(mean, std_dev)
        #plt.text(np.amax(n), mean+5, 'mean='+str(mean)+'\nstd_dev='+str(std_dev), horizontalalignment='left', fontsize=7)
        
        if np.amax(n) > n_max:
            n_max = np.amax(n)
    print(n_max)
    #ax_h.set_xlim(0, 22)
    #ax_h.set_xlim(0, n_max*4)
    ax_h.set_ylim(Ymin, Ymax)
    #ax_h.grid(True)
    plt.yticks([],[])
    ax_h.set_aspect(aspect = 0.3)
    plt.xticks([5, 10, 15, 20], ['5', '10', '15', '20'], fontsize=22)
    #plt.xlabel('number of cells')
    #plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', #'+str(Nrow)+'cells\n'+title, fontsize=7)
    plt.savefig(path_plot+plot_file+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_horizontal.pdf')
    plt.close()


def IDS_VGS(chip, col, L, Nfin, VT_flavor, Nrow, data_files, colors, path_plot, plot_file, row_idx, title, Ymax, IV_group_legend = [], VD = 0.8, VG_min = 0.2, VG_max = 0.8):

    if os.path.isdir(path_plot) == False:
        os.mkdir(path_plot)

    VDD_WL = np.arange(VG_max, VG_min -0.0001, -0.05)
    #Ymax = 0
#    for direction in ['VAsource_VBdrain', 'VAdrain_VBsource']:
    IV_group = []
    for (data_file, color) in zip(data_files, colors):
        Isense = []
        #f=open(path_data+'Fresh_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_Ids_Vgs_'+direction,'rU')
        f=open(data_file,'rU')
        Isense.append(re.findall(r'Isense=(-*\d*\.\d+)',f.read()))
        f.close()
        print(len(Isense), len(Isense[0]))
        if (len(Isense) != 1 or len(Isense[0]) != (Nrow+1)*len(VDD_WL)):
            print('data file error!')
        #I_leak_VAsource_VBdrain = np.zeros(len(VDD_WL))
        #I_leak_VAdrain_VBsource = np.zeros(len(VDD_WL))
        #for k in np.arange(len(VDD_WL)):
        #    I_leak_VAsource_VBdrain[k] = np.float64(Isense_VAsource_VBdrain[k][0])
        #    I_leak_VAdrain_VBsource[k] = np.float64(Isense_VAdrain_VBsource[k][0])
        IDS_VGS = np.zeros((Nrow, len(VDD_WL)))
        for row in np.arange(0, Nrow):
            IDS_VGS[row][:] = np.array(np.float64(Isense[0][(row+1)*len(VDD_WL): (row+2)*len(VDD_WL)]))

        #Ymax = max(1e6 * np.amax(IDS_VGS), Ymax)
        #plt.figure(figN)
        #plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', #'+str(Nrow)+' fresh nmos IDS-VGS curves', fontsize=10)
        for row in row_idx:
            #fig, = plt.plot(VDD_WL, 1e6*IDS_VGS[row], color=color, linestyle='solid', marker='.', alpha = 0.4, rasterized = True)
            fig, = plt.plot(VDD_WL, 1e6*IDS_VGS[row], color=color, linestyle='solid', marker='.', alpha = 0.4)
        IV_group.append(fig)

    #plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', '+title+' IDS-VGS curves', fontsize=7)
    if (len(IV_group_legend) != 0):
        plt.legend(IV_group, IV_group_legend, loc = 'best')
    plt.axis([VG_min, VG_max, 0, Ymax])
    plt.xlabel('$\mathregular{V_{GS}}$ (V)', fontsize = 23)
    plt.ylabel('$\mathregular{I_{DS}}$ ($\mathregular{\mu}$A)', fontsize = 23)
    #plt.grid()
    #plt.xlabel('VGS (V)')
    #plt.ylabel('IDS (uA)')
    plt.xticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], ['0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8'], fontsize=22)
    #plt.yticks([20, 40, 60, 80, 100], ['20', '40', '60', '80', '100'], fontsize=17)
    plt.yticks([20, 40, 60, 80, 100, 120, 140, 160], ['20', '40', '60', '80', '100', '120', '140', '160'], fontsize=22)
    #plt.grid()
    ##plt.draw()
    ##xticks = []
    ##xlocs, xlabels = plt.xticks()
    ##print(xlocs, xlabels)
    ##for xlabel in xlabels:
    ##    xticks.append(xlabel.get_text())
    ##    print(xlabel.get_text())
    ##print(xticks)
    ##yticks = []
    ##plt.xticks(xlocs, xticks, fontsize=17)
    ##ylocs, ylabels = plt.yticks()
    ##for ylabel in ylabels:
    ##    yticks.append(ylabel.get_text())
    ##plt.yticks(ylocs, yticks, fontsize=17)

    plt.show()
    plt.savefig(path_plot+plot_file+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'.pdf', dpi=300)
    #plt.savefig(path_plot+'IDS_VGS_'+str(figN).zfill(2)+'_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+direction+'.pdf')
    #figN = figN + 1
    plt.close()


if __name__ == '__main__':
  main()

