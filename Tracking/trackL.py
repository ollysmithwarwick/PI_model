import data
import os
import numpy as np

def outName(L):
    L = round(L*100)
    a = int(L/100)
    b = int((L-a*100)/10)
    c = int((L-a*100-b*10))
    name = 'Track' + str(a) + '_' + str(b) + str(c)
    return name

def numerateString(string):
    newString = ''.join(i for i in string if (i.isdigit()) or (i == '.'))
    return newString

def AddGridPoints():
    global gridpoints

    os.system('cp POPI_in0.dat POPI_in.dat')
    os.system(f'echo {gridpoints} | ./changegrid')
    os.system('cp POPI_in_newgrid.dat POPI_in0.dat')
    os.system('cp POPI_in1.dat POPI_in.dat')
    os.system(f'echo {gridpoints} | ./changegrid')
    os.system('cp POPI_in_newgrid.dat POPI_in1.dat')
    os.system('rm POPI_in_newgrid.dat')
    gridpoints += 2
    setup[3] = f' N_X = {gridpoints},'
    frameInfo[2] = f' boundR = {gridpoints},'
    data.WriteSetup(setup, 'setup.in')
    data.WriteFrameInfo(frameInfo, 'frameinfo.in')

def CleanUpGuess(guess):
    guess[1] = 0.0
    guess[2] = 0.0
    guess[4] = -1

def CopyTrackFiles():
    os.system(f"cp guess0.in ./{locName}")
    os.system(f"cp guess1.in ./{locName}")
    os.system(f"cp POPI_in0.dat ./{locName}")
    os.system(f"cp POPI_in1.dat ./{locName}")
    os.system(f"cp setup.in ./{locName}")

def CopyInputFiles():
    os.system(f"cp trackL.in ./{dirNameLong}")
    os.system(f"cp guess0_0.in ./{dirNameLong}")
    os.system(f"cp guess1_0.in ./{dirNameLong}")
    os.system(f"cp POPI_in0_0.dat ./{dirNameLong}")
    os.system(f"cp POPI_in1_0.dat ./{dirNameLong}")
    os.system(f'cp setup0.in ./{dirNameLong}')
    os.system(f'cp frameinfo0.in ./{dirNameLong}')
    os.system(f"cp trackPlotInfo.out ./{dirNameLong}")

def Run():
    os.system("../PIPO")

def Interpolate(step_, mode_):
    if (mode_ == 0):
        r = step_/(guessOut[0] - guessOutOld[0])
        nextGuess[0] = guessOut[0] + step_
    elif (mode_ == 1):
        r = 1
        nextGuess[0] = guessOut[0] + r * (guessOut[0] - guessOutOld[0])
        
    nextGuess[1] = 0.0
    nextGuess[2] = 0.0
    os.system(f'echo {r} | ./interp')
    
#Execute plan
## Loop over values of L consecutively
setup = []
setup0 = []
frameInfo = []
guessOut = []
guessOutOld = []
currentGuess = []
nextGuess = []
tempGuess = []
guess0 = []
guess1 = []
trackInfo = []

data.ReadGuess(guess0, 'guess0_0.in')
data.ReadGuess(guess1, 'guess1_0.in')
os.system('cp POPI_in0_0.dat POPI_in0.dat')
os.system('cp POPI_in1_0.dat POPI_in1.dat')
CleanUpGuess(guess0)
CleanUpGuess(guess1)
data.ReadSetup(setup, 'setup0.in')
data.ReadFrameInfo(frameInfo, 'frameinfo0.in')
os.system('cp frameinfo0.in frameinfo.in')
#data.ReadGuess(guess0, 'guess0.in')
#data.ReadGuess(guess1, 'guess1.in')
data.ReadTrackIn(trackInfo, 'trackL.in')

L0         = float(numerateString(setup[7]))
gridpoints =   int(numerateString(setup[3]))


mode = trackInfo[0]
stepS   = trackInfo[1]
n_S = trackInfo[2]
dirName = trackInfo[3]
dirNameLong = f'Tracking/{dirName}'
if not os.path.exists(dirNameLong):
    if not os.path.exists('./Tracking'):
        os.mkdir('./Tracking')
    os.mkdir(dirNameLong)


stepL = 0.05
dL = -stepL
n_L = 100
plotInfo = np.zeros(n_L)
for i_L in range(n_L):
    L = L0 + stepL*i_L
    dL = dL + stepL
    locName = f'{dirNameLong}/{outName(L)}'
    if not os.path.exists(f'{locName}'):
        os.mkdir(f'{locName}')
    if (round(dL*100) >= 50):
        AddGridPoints()
        dL = 0.0

    CleanUpGuess(guess0)
    CleanUpGuess(guess1)
    data.WriteGuess(guess0, 'guess0.in')
    data.WriteGuess(guess1, 'guess1.in')

    setup[1] = f" NAME = '{locName}/0',"
    setup[7] = f" L = {L}"

    data.WriteSetup(setup, 'setup.in')
    CopyTrackFiles()

    data.WriteGuess(guess0, 'guess.in')
    os.system('cp POPI_in0.dat POPI_in.dat')
    Run()
    guess0 = []
    data.ReadGuess(guess0, 'guess.out')
    guessOutOld = list(guess0)
    os.system("python plot2.py")
    os.system("cp out.dat out_old.dat")
    os.system("cp out.dat POPI_in0.dat")

    setup[1] = f" NAME = '{locName}/1',"
    data.WriteSetup(setup, 'setup.in')
    data.WriteGuess(guess1, 'guess.in')

    os.system('cp POPI_in1.dat POPI_in.dat')    
    Run()
    guess1 = []
    data.ReadGuess(guess1, 'guess.out')
    guessOut = list(guess1)
    os.system("cp out.dat POPI_in1.dat")
    os.system("python plot2.py")
    mode = 0
    nextGuess = list(guessOut)
    CleanUpGuess(nextGuess)

    for i_S in range(n_S):
        Interpolate(stepS, mode)
        os.system('cp out.dat out_temp.dat')
        guessTemp = list(guessOut)
        setup[1] = f" NAME = '{locName}/{i_S + 2 - mode}',"
        data.WriteSetup(setup, 'setup.in')
        data.WriteGuess(nextGuess, 'guess.in')
        Run()
        guessOut = []
        data.ReadGuess(guessOut, 'guess.out')
        CleanUpGuess(guessOut)
        if (guessOut[5] != 0):
            if (mode == 0):
                #            break
                mode = 1
                os.system('cp out_temp.dat out.dat')
                guessOut = list(guessTemp)
                continue
            else:
                print('mode 2 break')
                break
        else:
            guessOutOld = list(guessTemp)
            os.system('cp out_temp.dat out_old.dat')
            os.system("python plot2.py")
    
    plotInfo[i_L] = i_S

np.savetxt('trackPlotInfo.out', plotInfo)
CopyInputFiles()

### Take previous out.dat of runs 0 and 1 
### If L pushes dx too large, add 2 gridpoints
### 
#Print results
# Finish sorting out guess readin and out etc
# Finish sorting in and out dat files.
# Sort plotting properly
#Check Running
