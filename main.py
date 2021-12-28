import os
from tkinter import filedialog
from tkinter import *
import copy

global ACCEPTED
ACCEPTED = ['a','b','c','d','e','f','g']

global RESIDUES
RESIDUES = {
    "R":0,
    "H":0,
    "K":0,
    "D":0,
    "E":0,
    "S":0,
    "T":0,
    "N":0,
    "Q":0,
    "C":0,
    "G":0,
    "P":0,
    "A":0,
    "V":0,
    "I":0,
    "L":0,
    "M":0,
    "F":0,
    "Y":0,
    "W":0,
    "OTHER":0
}


class TALLEY:
    def __init__(self):
        self.count = 0
        self.flanks = [
            copy.deepcopy(RESIDUES),
            copy.deepcopy(RESIDUES),
            copy.deepcopy(RESIDUES)
        ]

    def Increment(self):
        self.count += 1


global FIRSTPOS
FIRSTPOS = {
    "R":TALLEY(),
    "H":TALLEY(),
    "K":TALLEY(),
    "D":TALLEY(),
    "E":TALLEY(),
    "S":TALLEY(),
    "T":TALLEY(),
    "N":TALLEY(),
    "Q":TALLEY(),
    "C":TALLEY(),
    "G":TALLEY(),
    "P":TALLEY(),
    "A":TALLEY(),
    "V":TALLEY(),
    "I":TALLEY(),
    "L":TALLEY(),
    "M":TALLEY(),
    "F":TALLEY(),
    "Y":TALLEY(),
    "W":TALLEY(),
    "OTHER":TALLEY()
}

def ScanTetrad(argv, sock, ind, val):
    Bools = [False, False, False]
    ind += 1
    found = 0
    missing = 0
    #Scan for argv[1] --> argv[0] --> argv[2]
    check = [1,0,2]
    group = list()
    group.append(val)
    try:
        while (sock.register[ind] != ' ' and sock.register[ind]) != '\n' or not False in Bools:
            if sock.register[ind] == argv[check[found]] and not Bools[argv.index(sock.register[ind])]:
                Bools[argv.index(sock.register[ind])] = True
                val = sock.sequence[ind]
                group.append(val)
                found = (found + 1) % 3
                missing = 0
            else:
                missing += 1
            if missing > 9:
                break
            ind += 1
    except IndexError:
        if not False in Bools:
            return group
    if not False in Bools:
        return group

def ScanTriad(argv, sock, ind, val):
    Bools = [False, True, False]
    ind += 1
    found = 0
    missing = 0
    #Scan for argv[0] --> argv[2]
    check = [0,2]
    group = list()
    group.append(val)
    try:
        while (sock.register[ind] != ' ' and sock.register[ind]) != '\n' or not False in Bools:
            if sock.register[ind] == argv[check[found]] and not Bools[argv.index(sock.register[ind])]:
                Bools[argv.index(sock.register[ind])] = True
                val = sock.sequence[ind]
                group.append(val)
                found = (found + 1) % len(check)
                missing = 0
            else:
                missing += 1
            if missing > 9:
                break
            ind += 1
    except IndexError:
        if not False in Bools:
            return group
    if not False in Bools:
        return group

def ScanTuple(argv, sock, ind, val):
    group = list()
    group.append(val)
    missing = 0
    while missing < 9 and ind < len(sock.register):
        if sock.register[ind] != argv[2]:
            missing += 1
            ind += 1
        else:
            group.append(sock.sequence[ind])
            return group

class SOCK_INFO:
    def __init__(self):
        self.register = None
        self.sequence = None
        self.flank = dict()

    def SetRegister(self, register):
        self.register = register[register.find(' '):] + " "

    def SetSequence(self, sequence):
        self.sequence = sequence[sequence.find(' '):] + " "

    def ContainsInfo(self):
        if self.register is not None and self.sequence is not None:
            return True
        return False


if __name__ == '__main__':
    print(sys.argv)
    argv = sys.argv[1:]
    if len(argv) != 3:
        sys.exit("Must have 3 parameters")
    for val in argv:
        if len(val) > 1 or val not in ACCEPTED:
            sys.exit(f"Invalid parameter: {val}")
    print(f"Checking flanking {argv[1]} and {argv[2]} residues using {argv[0]} as center")
    root = Tk()
    root.withdraw()
    folder_selected = filedialog.askdirectory()
    ind = 0
    tally = TALLEY()
    objects = list()
    tetrads = list()
    triads = list()
    tuples = list()
    for file in os.listdir(folder_selected):
        if file.endswith('.short.socket'):
            with open(folder_selected + '/'+ file) as sock_info:
                ind = ind +1
                pr = False
                currObj = SOCK_INFO()
                relLines = list()
                for line in sock_info:
                    if line.find('assigning heptad to helix') != -1:
                        pr = True
                    if line.isspace() or line.find('Finished') != -1:
                        pr = False
                        if currObj.ContainsInfo():
                            relLines.append(currObj)
                            currObj = SOCK_INFO()
                    if pr:
                        if line.find('sequence') == 0:
                            currObj.SetSequence(line)
                        if line.find('register') == 0:
                            currObj.SetRegister(line)
            #Retrieved all register-sequence pairs
            #Get flanking residues
            currCoil = dict()
            for structure in relLines:
                for i in range(len(structure.register)):
                    if structure.register[i] == argv[0]:
                        tetrad = ScanTetrad(argv, structure, i, structure.sequence[i])
                        tuple = ScanTuple(argv, structure, i, structure.sequence[i])
                        if tetrad is not None:
                            print("TETRAD (4): ", tetrad)
                            tetrads.append(tetrad)
                        if tuple is not None:
                            print("TUPLE  (2): ", tuple)
                            tuples.append(tuple)
                    elif structure.register[i] == argv[1]:
                        triad = ScanTriad(argv, structure, i, structure.sequence[i])
                        if triad is not None:
                            print("TRIAD  (3): ", triad)
                            triads.append(tetrad)

    print("Parsed everything")
    # outFile = open("data(old).csv", "w")
    # outFile.write(f"Data for {len(os.listdir(folder_selected))} socket files,\n")
    # outFile.write("\n")
    # for res in RESIDUES:
    #     outFile.write(f"\nFlanks for {res},\nFP Count,{FIRSTPOS[res].count},\n")
    #     outFile.write(",Residue:,")
    #     for res2 in RESIDUES:
    #         outFile.write(f"{res2},")
    #     outFile.write("\n,Next -b:,")
    #     for res2 in RESIDUES:
    #         outFile.write(f"{FIRSTPOS[res].flanks[0][res2]},")
    #     outFile.write(",\n")
    #     outFile.write(",Next -f:,")
    #     for res2 in RESIDUES:
    #         outFile.write(f"{FIRSTPOS[res].flanks[1][res2]},")
    #     outFile.write(",\n")
    #     outFile.write(",Next -c:,")
    #     for res2 in RESIDUES:
    #         outFile.write(f"{FIRSTPOS[res].flanks[2][res2]},")
    #     outFile.write("\n\n")
    # outFile.close()
