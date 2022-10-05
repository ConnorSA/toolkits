import numpy as np
import ase, ase.io
import re

class BandsFromCastep():
    RESCALE_TOL = 1e-5
    AU_TO_EV = 27.21132457027    
    def __init__(self, filein, kpath):
        self.read_in_file(filein)
        Klist_unsorted = self.get_Klist_unsorted()
        sort_idx, lines_diff =self.sort_Klist_get_idx(Klist_unsorted)
        eigenvec_unsorted = self.get_unsorted_eigenvecs(Klist_unsorted, lines_diff)
        self.sort_eigenvecs(sort_idx, eigenvec_unsorted)
        self.rescale_xaxis(kpath)
        return

    def read_in_file(self, filein):
        with open(filein) as f:
            self.filelines = f.readlines()
        return

    def get_fermi_level(self):
        for i, val in enumerate(self.filelines):
            if re.search('Fermi', val):
                break
        return float(self.filelines[i].split()[5])*self.AU_TO_EV

    def get_Klist_unsorted(self):
        return[ i for i, val in enumerate(self.filelines) if re.search(f'K-point', val) != None]
    def sort_Klist_get_idx(self, Klist_unsorted):
        if len(Klist_unsorted) < 2:
            return print('Only one k-point')
        lines_diff = Klist_unsorted[1] - Klist_unsorted[0]
        klist_raw = np.sort(self.filelines[Klist_unsorted[0]::lines_diff])
        self.kpath = np.array([i.split()[2:5] for i in klist_raw], dtype=float)
        return np.argsort(self.filelines[Klist_unsorted[0]::lines_diff]), lines_diff
    def get_unsorted_eigenvecs(self, Klist_unsorted, lines_diff):
        eigenvec_unsorted = []
        for i, val in enumerate(Klist_unsorted):
            eigenvec_unsorted.append([self.filelines[j] for j in range(val+2, val+lines_diff)])
        return np.array(eigenvec_unsorted, dtype=float)
    def sort_eigenvecs(self, sort_idx, eigenvec_unsorted):
        eigenvec_sorted=np.zeros_like(eigenvec_unsorted)
        for i, val in enumerate(sort_idx):
            eigenvec_sorted[i,:] = eigenvec_unsorted[val, :]
        self.eigenvec_sorted = eigenvec_sorted*AU_TO_EV
        return
    def find_index(self, kpath):
        j=0
        sympoint_idx=[]
        self.kpath_idx=[]
        for i, val in enumerate(self.kpath): #loop to find high symm points!
            if abs(np.linalg.norm(kpath[j]-val)) < self.RESCALE_TOL:
                sympoint_idx.append(i)
                j=j+1
        for i in range(len(sympoint_idx)-1):
            x = [j for j in range(sympoint_idx[i], sympoint_idx[i+1]+1)]
            self.kpath_idx.append(x)
    def rescale_xaxis(self, rescale_xaxis):
        #rescale axis to go 0 to 1 for equispaced kpath lines.
        in_path = np.array(rescale_xaxis, dtype=float)
        self.find_index(in_path)
        xsplit = 1/len(self.kpath_idx)
        xscale=[0]
        pos=0
        kpath_cut=[]
        kpath_cut.append(self.kpath_idx[0])
        for i in range(len(self.kpath_idx)-1):
            kpath_cut.append(self.kpath_idx[i+1][1:])
        x_inc = xsplit/(len(kpath_cut[0])-1)
        for j in range(len(kpath_cut[0])-1):
            pos = pos + x_inc
            xscale.append(pos)        
        for i in range(len(kpath_cut)-1):
            x_inc = xsplit/(len(kpath_cut[i+1]))
            for j in range(len(kpath_cut[i+1])):
                pos = pos + x_inc
                xscale.append(pos)
        self.xscale = np.array(xscale)




### example use hcp.
#note kpath should be identical to that specified in the castep.cell file.
#gamma = [0,0,0]
#K = [1/3,  1/3,  0]
#M = [0.5, 0, 0]
#A = [0, 0, 0.5]
#kpath = np.array([gamma, K, M, gamma, A])
#dz_0p01 = BandsFromCastep('location/to/file/castep.bands', kpath=kpath)

### Plotting
#labels = ["$\\Gamma$", "K", "M", "$\\Gamma$", "A"]
#ticks = np.linspace(0,1, len(labels))
#import matplotlib.pyplot as plt
#plt.plot(dz_0p01.xscale, dz_0p01.eigenvec_sorted)
#plt.hlines(xmin=0.05, xmax=1, y=dz_0p01.get_fermi_level(), color='k', linestyles='-')
#plt.text(0.01, dz_0p01.get_fermi_level(), r'$E_f$', ha='left', va='center', fontsize=16)
#plt.xticks(ticks, labels)
#plt.ylabel('Energy (eV)')
