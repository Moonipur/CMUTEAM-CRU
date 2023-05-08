#!/bin/python

import os, sys
from pathlib import Path
import pandas as pd
import numpy as np

def GEN_VAR():
    Input = sys.argv[1]
    Outdir = sys.argv[2]
    Outfile = sys.argv[3]

    return Input, Outdir, Outfile

def OUTDIR_CHECK(outdir):
    if os.path.isdir(f'{outdir}') == True:
        out = os.path.abspath(f'{outdir}')
    elif os.path.isdir(f'{outdir}') == False:
        print(IsADirectoryError)
        exit

    return out

def EXTRACT_CSV(inp):
    segment = pd.read_csv(f'{inp}')
    Chromosome = segment['Chromosome'].values.tolist()
    Start = segment['Start'].values.tolist()
    End = segment['End'].values.tolist()
    Copy = segment['Copies'].values.tolist()

    return Chromosome, Start, End, Copy

def ADD_CHR(Chromo):
    new_Chromo = []
    for i in Chromo:
        new = f'chr{int(i)}'
        new_Chromo.append(new)

    return new_Chromo

def DUP_or_DEL(copy):
    new_Type = []
    for i in copy:
        if i > 2:
            new = 'DUP'
        elif i < 2:
            new = 'DEL'
        elif i == 2:
            new = 'NOR'
        
        new_Type.append(new)
    
    return new_Type

def DELETE_NOR(n_Chr, start, end, n_Type):
    index_list = []
    for index, i in enumerate(n_Type, 0):
        if i == 'NOR':
            index_list.append(index)
        else:
            pass
    
    for j in sorted(index_list, reverse = True):
        del n_Chr[j]
        del start[j]
        del end[j]
        del n_Type[j]
    
    return n_Chr, start, end, n_Type

def MERGE(n_Chr, start, end, n_Type, outdir, outfile):
    merge_ = pd.DataFrame(
        list(zip(n_Chr, start, end, n_Type)),
        columns = ['CHR', 'START', 'END', 'TYPE']        
    )

    merge_.to_csv(f'{outdir}/{outfile}', index = False, header = False)

if __name__ == '__main__':
    VAR = GEN_VAR()
    OUT = OUTDIR_CHECK(VAR[1])
    DATA = EXTRACT_CSV(VAR[0])
    CHR = ADD_CHR(DATA[0])
    TYP = DUP_or_DEL(DATA[3])
    DEL = DELETE_NOR(CHR, DATA[1], DATA[2], TYP)
    MERGE(DEL[0], DEL[1], DEL[2], DEL[3], OUT, VAR[2])
