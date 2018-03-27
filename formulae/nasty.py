
from .utils import *

def nasty(mz,ppm=5):
    maxes = get_elemaxs(mz)
    error = mz *(ppm *1E-6)
    mlow, mhigh = mz-error, mz+error
    formula = {e:0 for e in eles.keys()}
    for H in range(maxes["H"]+1):
        formula["H"] = H
        mass = calc_mass(formula)
        if mass > mhigh:
            break
        if mass> mlow:
            yield formula.copy()
        for C in range(maxes["C"]+1):
            formula["C"] = C
            mass = calc_mass(formula)
            if mass > mhigh:
                break
            if mass> mlow:
                yield formula.copy()
            for N in range(maxes["N"]+1):
                formula["N"] = N
                mass = calc_mass(formula)
                if mass > mhigh:
                    break
                if mass> mlow:
                    yield formula.copy()
                for O in range(maxes["O"]+1):
                    formula["O"] = O
                    mass = calc_mass(formula)
                    if mass > mhigh:
                        break
                    if mass> mlow:
                        yield formula.copy()
                    for F in range(maxes["F"]+1):
                        formula["F"] = F
                        mass = calc_mass(formula)
                        if mass > mhigh:
                            break
                        if mass> mlow:
                            yield formula.copy()
                        for Si in range(maxes["Si"]+1):
                            formula["Si"] = Si
                            mass = calc_mass(formula)
                            if mass > mhigh:
                                break
                            if mass> mlow:
                                yield formula.copy()
                            for P in range(maxes["P"]+1):
                                formula["P"] = P
                                mass = calc_mass(formula)
                                if mass > mhigh:
                                    break
                                if mass> mlow:
                                    yield formula.copy()
                                for S in range(maxes["S"]+1):
                                    formula["S"] = S
                                    mass = calc_mass(formula)
                                    if mass > mhigh:
                                        break
                                    if mass> mlow:
                                        yield formula.copy()
                                    for Cl in range(maxes["Cl"]+1):
                                        formula["Cl"] = Cl
                                        mass = calc_mass(formula)
                                        if mass > mhigh:
                                            break
                                        if mass> mlow:
                                            yield formula.copy()
                                        for Br in range(maxes["Br"]+1):
                                            formula["Br"] = Br
                                            mass = calc_mass(formula)
                                            if mass > mhigh:
                                                break
                                            if mass> mlow:
                                                yield formula.copy()
                                            # elif mass >mlow:
                                            else:
                                                yield formula.copy()
