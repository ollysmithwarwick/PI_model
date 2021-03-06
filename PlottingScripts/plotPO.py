import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from shutil import copyfile

input = open("plotinfo.in")
name = (input.readline())
name = name.strip()
dirName = f'./PO/{name}/'
if not os.path.exists(dirName):
    os.mkdir(dirName)
    

copyfile('setup.in', f'{dirName}/setup.in')
#copyfile('init_mod.f90' , f'{dirName}/init.f90')
copyfile('frameinfo.in', f'{dirName}/frameinfo.in')
copyfile('init.dat' , f'{dirName}/init.dat' )
copyfile('POPI_in.dat', f'{dirName}/POPI_in.dat')
#copyfile('../../PI.f90', f'{dirName}/PI.f90')
copyfile('plotPO.py', f'{dirName}/plotPO.py')
#copyfile('../newton_PI.f90', f'{dirName}/newton_PI.f90')
copyfile('guess.in', f'{dirName}/guess.in')
copyfile('guess.out', f'{dirName}/guess.out')
copyfile('guesses.dat', f'{dirName}/guesses.dat')
copyfile('out.dat', f'{dirName}/out.dat')
copyfile('initFinal.dat', f'{dirName}/initFinal.dat')

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


flux = np.zeros((N_x,nframes+1))

xlist = np.linspace(0,nframes,nframes+1)
ylist = np.linspace(0,N_x-1,N_x)
x = np.arange(0,N_x,1)

X, Y = np.meshgrid(xlist,ylist)
#print(X)
#print(Y)

file = open("out.dat")
print(f'Printing Results for {nframes} frames')    
    
for t in np.arange(0, nframes+3):
#    print(t)

    nbar      = file.readline().strip().split()
    inbar     = file.readline().strip().split()
    E         = file.readline().strip().split()
    iE        = file.readline().strip().split()
    ntilde    = file.readline().strip().split()
    intilde   = file.readline().strip().split()
    phi       = file.readline().strip().split()
    iphi      = file.readline().strip().split()
    b_plus    = file.readline().strip().split()
    ib_plus   = file.readline().strip().split()
    b_minus   = file.readline().strip().split()
    ib_minus   = file.readline().strip().split()


    nbar      = np.array(nbar, 'float')
    E         = np.array(E, 'float')
    iE        = np.array(iE, 'float')
    inbar     = np.array(inbar, 'float')
    ntilde    = np.array(ntilde, 'float')
    intilde   = np.array(intilde, 'float')
    phi       = np.array(phi, 'float')
    iphi      = np.array(iphi, 'float')
    b_plus    = np.array(b_plus,'float')
    b_minus   = np.array(b_minus,'float')
    ib_plus   = np.array(ib_plus,'float')
    ib_minus  = np.array(ib_minus,'float')

#    for n in range(0,N_x):
#        flux[n,t] = (intilde[n]*phi[n]-iphi[n]*ntilde[n])
    if (t == 0):
        nbar0    = nbar
        E0       = E
        iE0      = iE
        inbar0   = inbar
        ntilde0  = ntilde
        intilde0 = intilde
        phi0     = phi
        iphi0    = iphi
        
    
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
        
   # ax5.plot(dx*x,b_plus[x])
   # ax5.plot(dx*x,ib_plus[x])
   # ax5.plot(dx*x, np.sqrt(b_plus[x]*b_plus[x]+ib_plus[x]*ib_plus[x]))
   # ax5.set_ylim(-2.0,2.0)
    
   # ax6.plot(dx*x, b_minus[x])
   # ax6.plot(dx*x, ib_minus[x])
   # ax6.plot(dx*x, np.sqrt(b_minus[x]*b_minus[x]+ib_minus[x]*ib_minus[x]))
   # ax6.set_ylim(-0.8,0.8)

#    ax6.plot(dx*x, flux[x,t])
#    ax6.set_ylim(-0.3,0.3)

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

    plt.tight_layout()
    if (t == nframes + 2):
        fig.savefig(f'./PO/{name}/final_shifted.png')
    elif (t == nframes + 1):
        fig.savefig(f'./PO/{name}/final.png')
    else:
        fig.savefig(f'./PO/{name}/fig{t}.png')
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
    

fig.savefig(f'./PO/{name}/residual.png')
plt.close()

#print(flux)
#fig2 = plt.figure()
#ax7 = fig2.add_subplot(111, projection = '3d')
#fig2, ax7 = plt.subplots(1)
#cp = ax7.contour(X, Y, flux)
#cp = ax7.plot_surface(X*spf*dt, Y*dx, flux)
#plt.colorbar(cp)
#fig2.savefig(f'./animation/{name}/contour.png')
#plt.show()    
#file.close()

