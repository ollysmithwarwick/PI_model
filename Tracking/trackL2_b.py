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

def AddGridPoints(gridpoints, setup, frameInfo):
    print(gridpoints)
    os.system('cp out.dat POPI_in.dat')
    os.system(f'echo {gridpoints} | ./changegrid')
    os.system('cp POPI_in_newgrid.dat out.dat')
    os.system('cp out_old.dat POPI_in.dat')
    os.system(f'echo {gridpoints} | ./changegrid')
    os.system('cp POPI_in_newgrid.dat out_old.dat')
    os.system('rm POPI_in_newgrid.dat')
    gridpoints += 2
    setup[3] = f' N_X = {gridpoints},'
    frameInfo[2] = f' boundR = {gridpoints},'
    data.WriteSetup(setup, 'setup.in')
    data.WriteFrameInfo(frameInfo, 'frameinfo.in')
    return gridpoints

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

def Interpolate(step_, mode_, guessOut, guessOutOld, nextGuess):
    if (mode_ == 0):
        r = step_/(guessOut[0] - guessOutOld[0])
        nextGuess[0] = guessOut[0] + step_
    elif (mode_ == 1):
        r = 1
        nextGuess[0] = guessOut[0] + r * (guessOut[0] - guessOutOld[0])
        
    nextGuess[1] = 0.0
    nextGuess[2] = 0.0
    os.system(f'echo {r} | ./interp')

def IsConverged(guess_):
    if(guess_[5] == 0):
        return True
    else:
        return False

# Read info in
# Run for initial L for n steps of length stepS
def main():
    setup      = []
    frameInfo  = []
    guess     = []
    nextGuess = []
    guessOut = []
    guessOutOld = []
    trackInfo  = []
    gridpoints = 0
    i_L = 0

    data.ReadSetup(setup, 'setup0.in')
    data.ReadFrameInfo(frameInfo, 'frameinfo0.in')
    data.ReadTrackIn(trackInfo, 'trackL.in')
    stepS = trackInfo[1]
    dirName = './Tracking/' + trackInfo[3]
    gridpoints0 = int(numerateString(setup[3]))
    for i_s in range(3):
        dirNameS = dirName + f'/{i_s}'
        gridpoints = gridpoints0
        setup[3] = f' N_X = {gridpoints},'
        frameInfo[2] = f' boundR = {gridpoints},'
        guess.clear()
        if i_s == 0:
            data.ReadGuess(guess, 'guess0.in')
            os.system('cp POPI_in0.dat POPI_in.dat')
        elif i_s == 1:
            data.ReadGuess(guess, 'guess1.in')
            os.system('cp POPI_in1.dat POPI_in.dat')
        else:
            os.system('cp out_oldL.dat out_old.dat')
            os.system('cp outL.dat out.dat')
            guessOut = list(guessOutL)
            guessOutOld = list(guessOutOldL)
            Interpolate(stepS, 0, guessOut, guessOutOld, nextGuess)
            guess = list(nextGuess)
            print(guessOut)
            print(guessOutOld)
            print(nextGuess)

        dL = -0.01

        for i_L in range(200):
            L = 5.0 + 0.01 * i_L
            dL = dL + 0.01
            if round(dL*100) >= 50:
                print(dL)
                print(gridpoints)
                gridpoints = AddGridPoints(gridpoints, setup, frameInfo)
                dL = 0.0
            if (i_L >= 2):
                Interpolate(1.0, 1, guessOut, guessOutOld, nextGuess)
                guessOutOld = list(guessOut)
                os.system('cp out.dat out_old.dat')
            elif (i_L == 0):
                nextGuess = list(guess)
            elif (i_L == 1):
                nextGuess = list(guessOut)
                guessOutOld = list(guessOut)
                os.system('cp out.dat POPI_in.dat')
                os.system('cp out.dat out_old.dat')

            dirNameL = str(dirNameS + '/' + outName(L))
            setup[1] = f" NAME ='{dirNameL}',"
            setup[7] = f' L = {L},'
            
            if not os.path.exists(dirNameL):
                if not os.path.exists(dirNameS):
                    if not os.path.exists(dirName):
                        os.mkdir(dirName)
                    os.mkdir(dirNameS)
                os.mkdir(dirNameL)

            print(nextGuess)
            CleanUpGuess(nextGuess)
            data.WriteGuess(nextGuess, 'guess.in')
            data.WriteSetup(setup, 'setup.in')
            data.WriteFrameInfo(frameInfo, 'frameinfo.in')
            Run()
            
            guessOut.clear()
            data.ReadGuess(guessOut, 'guess.out')
            if (IsConverged(guessOut)):
                os.system('python plot2.py')
                if i_L == 0:
                    if i_s != 0:
                        os.system('cp outL.dat out_oldL.dat')
                        guessOutOldL = list(guessOutL)
                    os.system('cp out.dat outL.dat')
                    guessOutL = list(guessOut)
                continue
            else:
                break
        if(i_L == 0):
            break

if (__name__ == '__main__'):
    main()
