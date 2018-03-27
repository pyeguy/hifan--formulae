from . import data
import os
from subprocess import run,PIPE
import io
import pandas as pd
import random
import re
from pyteomics import mass as pyteom_mass
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


cmd_path = os.path.abspath(os.path.join(__file__,"../../"))


def get_elemaxs(mz):
    top_masses = sorted(list(data.element_restrictions.keys()))
    for mass in top_masses:
        if mz <= mass:
            return data.element_restrictions[mass]
    return data.element_restrictions[top_masses[-1]]

simple_eles = {"C","H","O","N"}

def mass2formula_args(mz,ppm=5,all_ele=True,charge=None):

    tol = ((ppm * 1E-6 ) * mz) * (1E3)
    base_cmd = str(os.path.join(cmd_path,"formula.exe"))
    args_list = [base_cmd, f"-m {mz}",f"-t {tol}"]
    erestrict = get_elemaxs(mz)
    for ele,maxno in erestrict.items():
        if not all_ele:
            if ele in simple_eles:
                args_list.append(f"-{ele} 0-{maxno}")
            else:
                continue
        else:
            args_list.append(f"-{ele} 0-{maxno}")

    if charge:
        args_list.append(f"{charge}")
    return args_list

def run_formula(*args,**kwargs):
    args_list = mass2formula_args(*args,**kwargs)
    runo = run(args_list,stdout=PIPE)
    strout = runo.stdout.decode("utf-8")
    return strout

def formula_df(strout):
    sio = io.StringIO(strout)
    df = pd.read_table(sio,delimiter=';')
    return df

def fuzz_mass(mass,ppm=5,add_h=False):
    r = random.random()
    rppm = r*ppm
    rerror = rppm*1E-6*mass
    if r > 0.5:
        return mass + rerror
    else:
        return mass - rerror


def parse_formula(formula):
    part = re.findall(r'([A-Z][a-z]*)(\d*)',formula)
    done = [(atom,int(num)) if num else (atom,1) for atom,num in part]
    return done

def formula2emass(formula):
    em = pyteom_mass.calculate_mass(formula=formula)
    return em
    

megadb_path = r"C:\Users\camer\Desktop\free available databases in the review table 1\MegaDB\MegaDB.sdf"
def get_soome_mols(n=100):
    megadb = Chem.ForwardSDMolSupplier(megadb_path)
    some_mols = []
    for _ in range(100):
        m = next(megadb)
        if m:
            some_mols.append(m)
    return some_mols



eles = {
    "Br" : 78.9183371,
    "Cl" : 34.96885268,
    "S":  31.972071,
    "P" :  30.97376163,
    "Si" : 27.9769265325,
    "F" :  18.99840322,
    "O" :  15.9949146196,
    "N" :  14.0030740048,
    "C" :  12.000000000,
    "H" :   1.0078250321,
    }
elel = sorted([(s,m) for s,m in eles.items()],key = lambda x:x[1],reverse=True) #sort by mass
massl = [x[1] for x in elel]# only keep exactmass
elel = [x[0] for x in elel] # only keep symbol

def unparse_formula(form):
    strl = []
    for e,v in zip(elel,form):
        if v !=0:
            strl.append(f"{e}{v}")
    return ''.join(strl)


def calc_mass(form):
    mass = 0
    for m,cnt in zip(massl,form):
        if cnt != 0:
            mass += (m * cnt)
    return mass

def check_formula(form, fmass, mlow, mhigh):
    rval = False
    if mlow <= calc_mass(form) <= mhigh:
        rval = True
    return rval


def gen_formulae(mz,ppm=5,ratios=None):

    # elel = ["C","H","O"]
    # maxl = [3,3,3]
    maxes = dict(get_elemaxs(mz))
    maxl = [maxes[e] for e in elel]
    tmaxl = maxl.copy()
    error = mz *(ppm *1E-6)
    mlow, mhigh = mz-error, mz+error
    n = len(elel)
    formula = [0]*n
    i = 0
    bt = 0
    while formula[-1] < maxl[-1]:
        mass = calc_mass(formula)
        if check_formula(formula,mass,mlow,mhigh):
            yield formula.copy()
        i = 0 
        tmaxl = maxl.copy()
        while True:
            # print(formula.copy()," : ",i)
            if formula[i] < tmaxl[i]:
                formula[i] += 1
                if calc_mass(formula) > mhigh:
                    # do something here to backtrack...
                    if formula[i] <=0:
                        formula[i] = 0
                    else:
                        formula[i] -=1 
                    tmaxl[i] = formula[i]
                break

            else:
                formula[i] = 0
                i+=1
    

from itertools import product
def prod_soln(mz,ppm=5):
    maxes = dict(get_elemaxs(mz))
    maxl = [maxes[e] for e in elel]
    error = mz *(ppm *1E-6)
    mlow, mhigh = mz-error, mz+error
    formula = {e:0 for e in eles.keys()}
    pools = [[(ele,cnt) for cnt in range(maxes[ele]+1)] for ele in elel]
    pgen = product(*pools)
    for ft in pgen:
        formula = dict(ft)
        if check_formula(formula,mlow,mhigh):
            yield formula

def cartesian_product(aListOfList):
    indexes = [0] * len(aListOfList)
    while True:
        yield [l[i] for l,i in zip(aListOfList, indexes)]
        j = len(indexes) - 1
        while True:
            indexes[j] += 1
            if indexes[j] < len(aListOfList[j]): 
                break
            indexes[j] = 0
            j -= 1
            if j < 0: 
                return

def check_maxes(formd,maxes):
    bools = [v<maxes[e] for e,v in formd.items()]
    return all(bools)

def rec_formula(mz,ppm=5):
    maxes = dict(get_elemaxs(mz))
    error = mz *(ppm *1E-6)
    mlow, mhigh = mz-error, mz+error
    formula = {e:0 for e in eles.keys()}
    return _rec_form(formula,mlow,mhigh,maxes)

def _rec_form(ele_idx,formula,mlow,mhigh,maxes):
    good_form = check_formula(formula,mlow,mhigh)
    good_form = True
    under_maxes = check_maxes(formula,maxes)
    if good_form and under_maxes:
        yield formula
    else:
        formula[ele]
        pass






            
sm = get_soome_mols()
masses = [rdMolDescriptors.CalcExactMolWt(m) for m in sm]
formulas = [rdMolDescriptors.CalcMolFormula(m) for m in sm]

