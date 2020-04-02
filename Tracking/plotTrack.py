import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from shutil import copyfile
from shutil import copytree
from shutil import rmtree

input = open("Output/plotinfo.in")
name = (input.readline())
name = name.strip()
dirName = f'Tracking/Data/{name}'
if not os.path.exists(dirName):
    os.mkdir(dirName)
    
S = float(input.readline().strip())
nframes = (input.readline())
N_x = (input.readline())
L = (input.readline())
dt = input.readline()
dt = input.readline()
N_t = input.readline()
input.close()

nframes = nframes.strip()
N_x = N_x.strip()
L = L.strip()
N_t = N_t.strip()
dt = dt.strip()
nframes = int(nframes)
N_x = int(N_x)
L = float(L)
dt = float(dt)
N_t = int(N_t)

spf = int(N_t/nframes)
dx = L/N_x
total_t = N_t * dt
tpf = spf * dt

flux = np.zeros((N_x,nframes+1))

xlist = np.linspace(0,nframes,nframes+1)
ylist = np.linspace(0,N_x-1,N_x)
x = np.arange(0,N_x,1)

X, Y = np.meshgrid(xlist,ylist)
#print(X)
#print(Y)

file = open("POPI/Output/out.dat")
print(f'Printing Results for {nframes} frames')    
    
for t in np.arange(0, nframes+2):
#    print(t)
    time = t*tpf
    if(t >= nframes + 1):
        time = total_t

    nbar      = file.readline().strip().split()
    inbar     = file.readline().strip().split()
    E         = file.readline().strip().split()
    iE        = file.readline().strip().split()
    ntilde    = file.readline().strip().split()
    intilde   = file.readline().strip().split()
    phi       = file.readline().strip().split()
    iphi      = file.readline().strip().split()

    nbar      = np.array(nbar, 'float')
    E         = np.array(E, 'float')
    iE        = np.array(iE, 'float')
    inbar     = np.array(inbar, 'float')
    ntilde    = np.array(ntilde, 'float')
    intilde   = np.array(intilde, 'float')
    phi       = np.array(phi, 'float')
    iphi      = np.array(iphi, 'float')

    if (t <= nframes):
        for n in range(0,N_x):
            print(nframes, n, t)
            flux[n,t] = (intilde[n]*phi[n]-iphi[n]*ntilde[n])

    if (t == 0):
        nbar0    = nbar
        E0       = E
        iE0      = iE
        inbar0   = inbar
        ntilde0  = ntilde
        intilde0 = intilde
        phi0     = phi
        iphi0    = iphi
        
    if (t==0 or t == nframes + 2):
    # fig, ((ax1, ax2),(ax3,ax4),(ax5,ax6))  = plt.subplots(3,2)
        fig, ((ax1, ax2),(ax3,ax4))  = plt.subplots(2,2)
    
        ax1.plot(dx*x, phi[x])
        ax1.plot(dx*x, iphi[x])
        ax1.plot(dx*x, np.sqrt(phi[x]*phi[x]+iphi[x]*iphi[x]))
        ax1.set_ylim((-1.0,1.0))
    
        ax2.plot(dx*x,ntilde[x])
        ax2.plot(dx*x,intilde[x])
        ax2.plot(dx*x, np.sqrt(ntilde[x]*ntilde[x]+intilde[x]*intilde[x]))
        ax2.set_ylim(-1.0,1.0)
    
        ax3.plot(dx*x,nbar[x])
        ax3.plot(dx*x,inbar[x])
        ax3.set_ylim(-1.0,1.0)
        
        ax4.plot(dx*x,E[x])
        ax4.plot(dx*x,iE[x])
        ax4.set_ylim(-1.0,1.0)
        
        if (t == nframes and 1==0):
            ax1.set_title('$d\phi/dt$')
            ax2.set_title('$dE/dt$')
            ax3.set_title('$dn_0/dt$')
            ax4.set_title('$dn/dt$')
            #    ax5.set_title('$da_+/dt$')
            #    ax6.set_title('$da_-/dt$')
        else:
            ax1.set_title('$\phi$')
            ax2.set_title('$n$')
            ax3.set_title('$n_0$')
            ax4.set_title('$E$')
            #    ax5.set_title('$a_+$')
            #    ax6.set_title('$a_-$')
            #ax6.set_title('Flux')
        for a in [ax1, ax2, ax3, ax4]:
            a.set_xlabel('x')
            a.set_ylim(-1,1)

        fig.suptitle('t = ' + str(time)[0:7])
        plt.tight_layout()

        if (t == nframes + 2):
            fig.savefig(f'./Tracking/Data/{name}/final_shifted.png')
        elif (t == nframes + 1):
            fig.savefig(f'./Tracking/Data/{name}/final.png')
        else:
            fig.savefig(f'./Tracking/Data/{name}/fig{t}.png')
        plt.close()

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

ax1.plot(dx*x, phi[x]-phi0[x])
ax1.plot(dx*x, iphi[x]-iphi0[x])
ax2.plot(dx*x, ntilde[x]-ntilde0[x])
ax2.plot(dx*x, intilde[x]-intilde0[x])
ax3.plot(dx*x, nbar[x]- nbar0[x])
ax3.plot(dx*x, inbar[x]-inbar0[x])
ax4.plot(dx*x, E[x]-E0[x])
ax4.plot(dx*x, iE[x]-iE0[x])
ax1.set_title('Residual in $\phi$')
ax2.set_title('Residual in $n$')
ax3.set_title('Residual in $n0$')
ax4.set_title('Residual in $E$')

for a in [ax1, ax2, ax3, ax4]:
    a.set_xlabel('x')
plt.tight_layout()
    

fig.savefig(f'./{dirName}/residual.png')
plt.close()

np.savetxt(f'{dirName}/flux.dat', flux, fmt = '%10.6f')

peak_flux = np.amax(flux[:, :])
np.savetxt(f'{dirName}/peakflux.dat', np.array([S, peak_flux]), fmt = '%15.8f')


#if os.path.exists(f'{dirName}/Input'):
#    rmtree(f'{dirName}/Input')
#copytree('Input', f'{dirName}/Input')
copyfile('PlottingScripts/plot.py', f'{dirName}/plot.py')
#if os.path.exists(f'{dirName}/Output'):
#    rmtree(f'{dirName}/Output')
#copytree('Tracking/Output/', f'{dirName}/Output/')


fig2, ax2 = plt.subplots()
cs = ax2.contour(X, Y, flux, 20)
plt.tight_layout()
fig2.savefig(f'./{dirName}/contour.png')


