import fileinput

def ReadTrackIn(listin, filename):
    filein = open(filename, 'r')
    listin.append(  int(filein.readline().strip('\n')))
    listin.append(float(filein.readline().strip('\n')))
    listin.append(  int(filein.readline().strip('\n')))
    listin.append(  str(filein.readline().strip('\n')))
    filein.close()
    
def ReadSetup(listin, filename): # Since it is a namelist each line stored as string - have to manually manipulate
    filein = open(filename, 'r')

    listin.append(str(filein.readline().strip('\n')))
    listin.append(str(filein.readline().strip('\n')))
    listin.append(str(filein.readline().strip('\n')))
    listin.append(str(filein.readline().strip('\n')))
    listin.append(str(filein.readline().strip('\n')))
    listin.append(str(filein.readline().strip('\n')))
    listin.append(str(filein.readline().strip('\n')))
    listin.append(str(filein.readline().strip('\n')))
    listin.append(str(filein.readline().strip('\n')))
    listin.append(str(filein.readline().strip('\n')))

    filein.close()

def ReadFrameInfo(listin, filename):
    
    filein = open(filename, 'r')
    for i in range(8):
        listin.append(str(filein.readline().strip('\n')))

    filein.close()

def ReadGuess(listin, filename):
    filein = open(filename, 'r')
    listin.append(float(filein.readline().strip('\n'))) # S
    listin.append(float(filein.readline().strip('\n'))) # L shift
    listin.append(float(filein.readline().strip('\n'))) # phase shift
    listin.append(  int(filein.readline().strip('\n'))) # N_BCperiods
    listin.append(  int(filein.readline().strip('\n'))) # N_timesteps
    listin.append(  int(filein.readline().strip('\n'))) # Info
    listin.append(  int(filein.readline().strip('\n'))) # N_newtonsteps
    filein.close()

def WriteSetup(listin, filename):
    fileout = open(filename, 'w')
    for i in range(10):
        a = listin[i]
        fileout.write(f'{a}\n')

    fileout.close()

def WriteTrackIn(listin, filename):
    fileout = open(filename, 'w')
    for i in range(4):
        a = str(listin[i])
        fileout.write(f"{a}\n")

    fileout.close()
    
def WriteFrameInfo(listin, filename):
    fileout = open(filename, 'w')
    for i in range(8):
        a = listin[i]
        fileout.write(f'{a}\n')
        
    fileout.close()

#def WriteFrameInfo(listin, filename):

def WriteGuess(listin, filename):
    fileout = open(filename, 'w')
    for i in range(3):           # S, L_shift, phase_shift
        a = float(listin[i])
        fileout.write(f"{a}\n")

    for i in range(3,7):         # N_BC_periods, N_timesteps, Info, N_newtonsteps
        a = int(listin[i])
        fileout.write(f"{a}\n")

    fileout.close()
    
def main():
    setup = []
    ReadSetup(setup, 'setup.in')
    

if __name__ == '__main__':
    main()
