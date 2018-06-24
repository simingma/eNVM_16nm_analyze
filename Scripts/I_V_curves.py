
import sys
import re
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import interp1d

def hist_IDS_VGS_baking(VG_idx, chip, col, L, Nfin, VT_flavor, Nrow, baking_data_files, bake_times, colors, legend_levels, path_plot, plot_file, row_idx, title, Ymin, Ymax, Imin, Imax, Num_bins = 50, VD = 0.8, VG_min = 0.2, VG_max = 0.8):
    """modified from hist_IDS_VGS to account for changes during baking"""

    if os.path.isdir(path_plot) == False:
        os.mkdir(path_plot)

    VDD_WL = np.arange(VG_max, VG_min -0.0001, -0.05)
    #Ymax = 0
#    for direction in ['VAsource_VBdrain', 'VAdrain_VBsource']:
    fig_h, ax_h = plt.subplots(nrows = 1, ncols = 1)
    n_max = 0
    #mean = np.zeros(len(data_files))
    #std_dev = np.zeros(len(data_files))
    mean_bake = np.zeros((len(baking_data_files), len(baking_data_files[0])))
    std_bake = np.zeros((len(baking_data_files), len(baking_data_files[0])))
    for (time_idx, data_files) in zip(np.arange(len(baking_data_files)), baking_data_files):
        for (level_idx, data_file, color) in zip(np.arange(len(data_files)), data_files, colors):
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
            mean_bake[time_idx][level_idx] = mean
            std_bake[time_idx][level_idx] = std_dev
            #plt.text(np.amax(n), mean+5, 'mean='+str(mean)+'\nstd_dev='+str(std_dev), horizontalalignment='left', fontsize=7)
            
            if np.amax(n) > n_max:
                n_max = np.amax(n)
        #print(n_max)
        #ax_h.set_xlim(0, 40)
        ##ax_h.set_xlim(0, n_max*4)
        #ax_h.set_ylim(Ymin, Ymax)
        ##ax_h.grid(True)
        #plt.yticks([],[])
        ##plt.xlabel('number of cells')
        ##plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', #'+str(Nrow)+'cells\n'+title, fontsize=7)
        #plt.savefig(path_plot+plot_file+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_horizontal.pdf')
        #plt.close()

        """
        fig_v, ax_v = plt.subplots(nrows = 1, ncols = 1)
        n_max = 0
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
            plt.figure(2)
            n, bins, pathes = ax_v.hist(1e6*IDS_VGS[row_idx, [VG_idx]*len(row_idx)], bins=Num_bins, normed = False, range=(Imin, Imax), orientation = 'vertical', color = color, edgecolor='none')
            if np.amax(n) > n_max:
                n_max = np.amax(n)
        print(n_max)
        #ax_v.set_ylim(0, n_max)
        ax_v.grid(True)
        plt.ylabel('number of cells')
        plt.xlabel('IDSAT (uA)')
        ax_v.set_xlim(Ymin, Ymax)
        ax_v.set_ylim(0, n_max)
        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', #'+str(Nrow)+'cells\n'+title, fontsize=7)
        plt.savefig(path_plot+plot_file+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_vertical.pdf')
        plt.close()
        """
    plt.close()
    plt.subplot(121)
    legend = []
    for (color, level_idx) in zip(colors, np.arange(len(baking_data_files[0]))):
        axe, = plt.plot(bake_times, mean_bake[:, level_idx], color=color, marker = '.')
        print('level='+str(level_idx)+'\n')
        for I_t in mean_bake[:, level_idx]:
            print((I_t-mean_bake[0, level_idx])/mean_bake[0, level_idx])
        legend.append(axe)
    plt.ylim(ymax = 120)
    #plt.legend(legend, legend_levels)
    plt.xticks([0, 20, 40, 70], ['0', '20', '40', '70'], fontsize=17)
    plt.yticks([0, 20, 40, 60, 80, 100, 120], ['0', '20', '40', '60', '80', '100', '120'], fontsize=17)
    #axe.set_aspect(aspect=0.5)
    #plt.savefig(path_plot+plot_file+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_mean_baking.pdf')
    #plt.close()
    plt.subplot(122)
    for (color, level_idx) in zip(colors, np.arange(len(baking_data_files[0]))):
        axe, = plt.plot(bake_times, std_bake[:, level_idx], color=color, marker = '.')
    plt.ylim(ymin=0)
    plt.xticks([0, 20, 40, 70], ['0', '20', '40', '70'], fontsize=17)
    plt.yticks([0, 0.4, 0.8, 1.2, 1.6, 2.0], ['0', '0.4', '0.8', '1.2', '1.6', '2.0'], fontsize=17)
    #axe.set_aspect(aspect=0.5)
    #plt.savefig(path_plot+plot_file+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_std_baking.pdf')
    plt.show()
    plt.savefig(path_plot+plot_file+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_mean__std_baking.pdf')
    plt.close()


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

    """
    fig_v, ax_v = plt.subplots(nrows = 1, ncols = 1)
    n_max = 0
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
        plt.figure(2)
        n, bins, pathes = ax_v.hist(1e6*IDS_VGS[row_idx, [VG_idx]*len(row_idx)], bins=Num_bins, normed = False, range=(Imin, Imax), orientation = 'vertical', color = color, edgecolor='none')
        if np.amax(n) > n_max:
            n_max = np.amax(n)
    print(n_max)
    #ax_v.set_ylim(0, n_max)
    ax_v.grid(True)
    plt.ylabel('number of cells')
    plt.xlabel('IDSAT (uA)')
    ax_v.set_xlim(Ymin, Ymax)
    ax_v.set_ylim(0, n_max)
    plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', #'+str(Nrow)+'cells\n'+title, fontsize=7)
    plt.savefig(path_plot+plot_file+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_vertical.pdf')
    plt.close()
    """

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

def IDS_VGS_stress(chip, col, L, Nfin, VT_flavor, Nrow, VGS, VDS, Nhci, Npbti, VD = 0.8, VG_min = 0.2, VG_max = 0.8, PulseCycle = 3, stress_time=['400ms','800ms','1200ms']):

    path_data = '../data/VG_ConstPulse_chip'+str(chip).zfill(2)+'/'
    path_plot = '../plot/VG_ConstPulse_chip'+str(chip).zfill(2)+'_IDS_VGS/'

    if os.path.isdir(path_plot) == False:
        os.mkdir(path_plot)

    VDD_WL = np.arange(VG_max, VG_min -0.0001, -0.05)
    for direction in ['VAsource_VBdrain', 'VAdrain_VBsource']:
        Isense = []
        f=open(path_data+'Fresh_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_Ids_Vgs_'+direction,'rU')
        Isense.append(re.findall(r'Isense=(-*\d*\.\d+)',f.read()))
        f.close()
        print(len(Isense), len(Isense[0]))
        if (len(Isense) != 1 or len(Isense[0]) != (Nrow+1)*len(VDD_WL)):
            print('data file error!')
        #I_leak_VAsource_VBdrain = np.zeros(len(VDD_WL))
        #for k in np.arange(len(VDD_WL)):
        #    I_leak_VAsource_VBdrain[k] = np.float64(Isense_VAsource_VBdrain[k][0])
        IDS_VGS = np.zeros((Nrow, len(VDD_WL)))
        for row in np.arange(0, Nrow):
            IDS_VGS[row][:] = np.array(np.float64(Isense[0][(row+1)*len(VDD_WL): (row+2)*len(VDD_WL)]))
        
        IDS_VGS_stress = np.zeros((PulseCycle, Nrow, len(VDD_WL)))
        for cycle in np.arange(PulseCycle):
            Isense = []
            f=open(path_data+'Stress_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_Ids_Vgs_'+direction+'_'+str(cycle+1).zfill(2),'rU')
            Isense.append(re.findall(r'Isense=(-*\d*\.\d+)',f.read()))
            f.close()
            print(len(Isense), len(Isense[0]))
            if (len(Isense) != 1 or len(Isense[0]) != (Nrow+1)*len(VDD_WL)):
                print('data file error!')
            for row in np.arange(0, Nrow):
                IDS_VGS_stress[cycle,row,:] = np.array(np.float64(Isense[0][(row+1)*len(VDD_WL): (row+2)*len(VDD_WL)]))


        Ymax = 1e6 * np.amax([np.amax(IDS_VGS), np.amax(IDS_VGS_stress)])
        #plt.figure(figN)
        plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', fresh nmos IDS-VGS curves ('+direction+')', fontsize=10)
        row = 0
        while(row < Nrow):
            for hci in np.arange(0, Nhci):
                fresh, = plt.plot(VDD_WL, 1e6*IDS_VGS[row], color='b', linestyle='solid', marker='.')
                row += 1
            for pbti in np.arange(0, Npbti):
                row += 1

        plt.axis([VG_min, VG_max, 0, Ymax])
        plt.xlabel('VGS (V)')
        plt.ylabel('IDS (uA)')
        plt.grid()
        plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+plot_file+'.pdf')
        plt.close()

        for cycle in np.arange(PulseCycle):
            plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', IDS-VGS curves ('+direction+') of fresh and '+stress_time[cycle]+' HCI\n VDS='+str(VDS)+'V, VGS='+str(VGS)+'V', fontsize=10)
            row = 0
            while(row < Nrow):
                for hci in np.arange(0, Nhci):
                    fresh, = plt.plot(VDD_WL, 1e6*IDS_VGS[row], color='b', linestyle='solid', marker='.')
                    stressed, = plt.plot(VDD_WL, 1e6*IDS_VGS_stress[cycle,row], color='r', linestyle='solid', marker='.')
                    row += 1
                for pbti in np.arange(0, Npbti):
                    row += 1

            plt.legend([fresh, stressed], ['fresh', 'HCI stress'], loc = 'best')
            plt.axis([VG_min, VG_max, 0, Ymax])
            plt.xlabel('VGS (V)')
            plt.ylabel('IDS (uA)')
            plt.grid()
            plt.savefig(path_plot+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+direction+'_'+str(cycle+1)+'.pdf')
            plt.close()


def IDSAT_vs_row(chip, col, L, Nfin, VT_flavor, Nrow, data_files, colors, path_plot, plot_file, row_idx, title, Ymax, IV_group_legend = [], VD = 0.8, VG_min = 0.2, VG_max = 0.8):

    """ Directly adapted from IDS_VGS, plot IDSAT vs row index (location) to check whether there's any systematic variation (indicative of IR drop), and also check the cell-to-cell variation"""

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
        fig, = plt.plot(row_idx, 1e6*IDS_VGS[row_idx,[0]*len(row_idx)], color=color, linestyle='solid', marker='.', alpha = 1.0)
        IV_group.append(fig)

    plt.title('L='+str(L)+', Nfin='+str(Nfin)+', '+VT_flavor+', '+title+' IDSAT vs row index/location', fontsize=7)
    if (len(IV_group_legend) != 0):
        plt.legend(IV_group, IV_group_legend, loc = 'best', fontsize = 7)
    #plt.axis([VG_min, VG_max, 0, Ymax])
    plt.xlim(min(row_idx), max(row_idx))
    plt.xlabel('row index')
    plt.ylabel('IDSAT (uA)')
    plt.grid()
    plt.savefig(path_plot+plot_file+'Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'.pdf')
    #plt.savefig(path_plot+'IDS_VGS_'+str(figN).zfill(2)+'_Chip'+str(chip).zfill(2)+'_Col'+str(col).zfill(2)+'_'+direction+'.pdf')
    plt.close()


if __name__ == '__main__':
  main()

