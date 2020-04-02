import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import trackL2
from shutil import copyfile
import os

def main():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    for iS in range(1, 200):
        S = []
        Amps = []
        L_list    = []
        usedSolns = []
        for iL in range(200):
            L = 5.0 - iL*0.05
            if (os.path.exists(f'./Tracking/Data/test_neg/{iS}/{trackL2.outName(L)}/peakflux.dat')):
                usedSolns.append(f'{trackL2.outName(L)}')
                L_list.append(L)
            else:
                break
        
        L_list.reverse()
        usedSolns.reverse()
        read(usedSolns, S, Amps, f'./Tracking/Data/test_neg/{iS}')
    
        usedSolns = []
        for iL in range(200):
            L = 5.0 + iL*0.02
            if (os.path.exists(f'./Tracking/Data/test_highS3/{iS}/{trackL2.outName(L)}/peakflux.dat')):
                usedSolns.append(f'{trackL2.outName(L)}')
                L = L_list.append(L)
            else:
                break
        read(usedSolns, S, Amps, f'./Tracking/Data/test_highS3/{iS}')
        usedSolns = []
        for iL in range(200):
            L = 5.0 + iL*0.05
            if (os.path.exists(f'./Tracking/Data/test6/{iS}/{trackL2.outName(L)}/peakflux.dat')):
                usedSolns.append(f'{trackL2.outName(L)}')
                L = L_list.append(L)
            else:
                break
        read(usedSolns, S, Amps, f'./Tracking/Data/test6/{iS}')

        usedSolns = []
        for iL in range(200):
            L = 5.0 + iL*0.05
            if (os.path.exists(f'./Tracking/Data/test_lowS3/{iS}/{trackL2.outName(L)}/peakflux.dat')):
                usedSolns.append(f'{trackL2.outName(L)}')
                L = L_list.append(L)
            else:
                break
        read(usedSolns, S, Amps, f'./Tracking/Data/test_lowS3/{iS}')

        plot(ax, S, L_list, Amps, f'L = {L}, Main Branch')
        
#    for i in range(10):
#        L = 5.0 + i*0.05
#        usedSolns.append(f'{trackL2.outName(L)}')

#    print(usedSolns)
#    read(usedSolns, S, Amps, 2)
#    plot(S, Amps, 'L = 5.0, Main Branch')

    ax.set_xlabel('S')
    ax.set_ylabel('Length of Box')
    ax.set_zlabel('Peak Amplitude')
#    plt.legend()
    plt.show()

def read(names, out1, out2, direc):
    for name in names:
        reads = []
        if os.path.exists(f'{direc}/{name}/peakflux.dat'):
            file = open(f'{direc}/{name}/peakflux.dat')
            reads = file.read().split()
            file.close()

            out1.append(float(reads[0]))
            out2.append(float(reads[1]))
            
def plot(ax, x, y, z, label_):
#    fig = plt.figure()
#    plt.plot(x,y, marker = 'x', label = label_)
#    x = np.array(x)
#    y = np.array(y)
#    z = np.array(z)
    ax.scatter(x,y, z, marker = 'x', label = label_)
    #    ax.plot_surface(x,y,z, cmap = 'viridis', edgecolor = 'none')
    x.clear()
    y.clear()

def copyFiles(dirName):
    copyfile('track2.in', f'{dirName}/track2.in')
    

if __name__ == '__main__':
    main()
