import numpy as np
import matplotlib.pyplot as plt
import os
from shutil import copyfile

input = open("plotinfo.in")
name = (input.readline())
S = (input.readline())
nframes = (input.readline())
N_x = (input.readline())
L = (input.readline())
time_per_frame = (input.readline())
dt = input.readline()
N_t = (input.readline())
input.close()

name = name.strip()
S = S.strip()
nframes = nframes.strip()
N_x = N_x.strip()
L = L.strip()
time_per_frame = float(time_per_frame.strip())
nframes = int(nframes)
N_x = int(N_x)
L = float(L)
dt = float(dt.strip())
N_t = int(N_t.strip())
dx = L/N_x
a = 0

fig, (ax1)  = plt.subplots(1,1)

file = open("bis.dat")
for line in file:
    a = a+1
    amp = (line.strip().split())
    amp = [float(i) for i in amp]
#    print(a, len(amp))
#    print(amp)
    while (amp[-1] == 0):
        amp.remove(0)

#    print(len(amp))
    
    t = np.arange(0, len(amp))
    if (a <= 54):
        ax1.semilogy(time_per_frame*t, amp, c = 'red',linewidth = 0.5)
    elif (a<=55):
        ax1.semilogy(time_per_frame*t, amp, c = 'green')
    else:
        ax1.semilogy(time_per_frame*t, amp, c = 'blue')

ax1.set_xlabel('t')
#ax1.set_xlim(-10,180)
ax1.set_ylabel('Amplitude')
#ax1.set_ylim(10^(-3),1.0)
plt.tight_layout()
#plt.show()
if( not(os.path.isdir(f'bisection'))):
    os.mkdir(f'bisection')
if( not(os.path.isdir(f'bisection/{name}'))):
    os.mkdir(f'bisection/{name}')

copyfile('bis.dat', f'bisection/{name}/bis.dat')
copyfile('bisect.in', f'bisection/{name}/bisect.in')
copyfile('setup.in', f'bisection/{name}/setup.in')
copyfile('bisect.f90', f'bisection/{name}/bisect.f90')
copyfile('PI.f90', f'bisection/{name}/PI.f90')
copyfile('out.dat', f'bisection/{name}/out.dat')
copyfile('plotBis.py', f'bisection/{name}/plotBis.py')
fig.savefig(f'./bisection/{name}/amp.png')

plt.close()

file.close()

file = open("out.dat")
print('READ')
flux = np.zeros((N_x, nframes + 1))
    
for t in np.arange(0, nframes+1):

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

    for n in range(0,N_x):
        flux[n,t] = (intilde[n]*phi[n]-iphi[n]*ntilde[n])
    
    x = np.arange(0,N_x,1)
    
#    fig, ((ax1, ax2),(ax3,ax4),(ax5,ax6))  = plt.subplots(3,2)
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
        
    #ax5.plot(dx*x,b_plus[x])
    #ax5.plot(dx*x,ib_plus[x])
    #ax5.plot(dx*x, np.sqrt(b_plus[x]*b_plus[x]+ib_plus[x]*ib_plus[x]))
    #ax5.set_ylim(-2.0,2.0)
    
    # ax6.plot(dx*x, b_minus[x])
    # ax6.plot(dx*x, ib_minus[x])
    # ax6.plot(dx*x, np.sqrt(b_minus[x]*b_minus[x]+ib_minus[x]*ib_minus[x]))
    # ax6.set_ylim(-0.8,0.8)

    #ax6.plot(dx*x, flux[x,t])
    #ax6.set_ylim(-0.3,0.3)

    if (t == nframes and 1==0):
        ax1.set_title('$d\phi/dt$')
        ax2.set_title('$dE/dt$')
        ax3.set_title('$dn_0/dt$')
        ax4.set_title('$dn/dt$')
       # ax5.set_title('$da_+/dt$')
       # ax6.set_title('$da_-/dt$')
    else:
        ax1.set_title('$\phi$')
        ax2.set_title('$n$')
        ax3.set_title('$n_0$')
        ax4.set_title('$E$')
       # ax5.set_title('$a_+$')
       # ax6.set_title('$a_-$')
       # ax6.set_title('Flux')

    plt.tight_layout()
    fig.savefig(f'./bisection/{name}/fig{t}.png')
    plt.close()
