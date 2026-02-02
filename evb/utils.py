import copy
import openmm as mm
import openmm.app as app
import numpy as np
from openmm.unit import *
from numba import *


def create_Mysystem(psf,ff,box_wu,nbcut,switchdist,ewaldtol,isPBC,isNVE,isoff):
    if isPBC:
        psf.setBox(box_wu[0],box_wu[1],box_wu[2])
        system = ff.createSystem(psf.topology,nonbondedMethod=app.PME,nonbondedCutoff=nbcut,switchDistance=switchdist,ewaldErrorTolerance=ewaldtol,constraints=None,rigidWater=False)
        modify_nbfix(system, nbcut, switchdist)
        for i in range(system.getNumForces()):
        # #     # myforce = system.getForce(i)
        # #     # forceName=myforce.__class__.__name__
        # #     # print("Force name: %s, Periodic: %s"%(forceName,myforce.usesPeriodicBoundaryConditions()))
        # #     # if "Harmonic" in forceName or "Torsion" in forceName:
        # #     #     print("Converting this Force %s to use PBC distances"%forceName)    
        # #     #     myforce.setUsesPeriodicBoundaryConditions(True)
            myforce = system.getForce(i)
            # if not myforce.usesPeriodicBoundaryConditions():
            #     try:
            #         print("##Setting %s PBC to True for this force"%(myforce.__class__.__name__))
            #         myforce.setUsesPeriodicBoundaryConditions(True)
            #     except:
            #         pass
                
            # if not myforce.usesPeriodicBoundaryConditions():
            #     print("##WARNING: can't set PBC to True for this force")
            #     print(myforce.__class__.__name__)
                    
    else:
        system = ff.createSystem(psf.topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    if isNVE and not isoff: add_drude(system, psf)
    return system

def wrappos(teststate,box,atomlist,targetindex):
    testposL=teststate.getPositions(asNumpy=True)._value            
    targetpos=testposL[targetindex]
    for index in atomlist:
       testposL[index]-=np.round((testposL[index]-targetpos)/box)*box 
    return testposL
           

def make_simulation(psf, system, Temp, dt,isNVE):
    platform = mm.Platform.getPlatformByName('CUDA')
    platform.loadPluginsFromDirectory('/home/yutong/software/myOpenMM-plumed-plugin/openmm-plumed1.0_install/lib/plugins/')
    if not isNVE: ## it is NVT
        integrator = mm.DrudeLangevinIntegrator(Temp*kelvin, 5/picosecond, 1*kelvin, 20/picosecond, dt)
        integrator.setMaxDrudeDistance(0.02)
        integrator.setConstraintTolerance(1e-06)
        prop = dict(CudaPrecision='mixed')
    else: ## it is NVE
        integrator = mm.DrudeSCFIntegrator(dt)
        integrator.setMinimizationErrorTolerance(1e-05)
        prop = dict(CudaPrecision='double')

    simulation = app.Simulation(psf.topology, system, integrator, platform, prop)
    return simulation


def add_drude(system, psf):
    hyper = mm.CustomBondForce('(step(r-rhyper)*(r-rhyper)*khyper)^powh')
    hyper.addGlobalParameter('khyper', 100.0)
    hyper.addGlobalParameter('rhyper', 0.020)
    hyper.addGlobalParameter('powh', 6)
    system.addForce(hyper)
    for d in psf.drudepair_list:
        hyper.addBond(d[0], d[1], [])

def modify_nbfix(system, nbcut, switchdist):
    nbfix = [
        system.getForce(i) for i in range(system.getNumForces())
        if isinstance(system.getForce(i), mm.CustomNonbondedForce)
    ][0]
    nbfix.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    nbfix.setCutoffDistance(nbcut)
    nbfix.setUseSwitchingFunction(True)
    nbfix.setSwitchingDistance(switchdist)

def check_force(system):
    nbforce = [
        system.getForce(i) for i in range(system.getNumForces())
        if isinstance(system.getForce(i), mm.NonbondedForce)
    ][0]

    print(nbforce.getNumParticles())

    print(nbforce.getSwitchingDistance())
    print(nbforce.getCutoffDistance())
    print(nbforce.usesPeriodicBoundaryConditions())

    for i in range(10):
        print(nbforce.getParticleParameters(i))

    nbfix = [
        system.getForce(i) for i in range(system.getNumForces())
        if isinstance(system.getForce(i), mm.CustomNonbondedForce)
    ][0]

    print(nbfix.getNumParticles())
    print(nbfix.getSwitchingDistance())
    print(nbfix.getCutoffDistance())
    print(nbfix.usesPeriodicBoundaryConditions())

    for i in range(10):
        print(nbfix.getParticleParameters(i))

def forcegroupify(system):
  #  wantedForceNames=["HarmonicBondForce","NonbondedForce","DrudeForce","HarmonicAngleForce","CustomNonbondedForce","CustomCVForce",]
    forcegroups = {}
    for i in range(system.getNumForces()):
        force = system.getForce(i)     
  #      if force.__class__.__name__  in wantedForceNames:
        force.setForceGroup(i)
        forcegroups[force] = i
    return forcegroups

def printDecomposedEnergies(state, simulation, fgrps,printHeaders):    
    forcenames="PE    "               
    values="%.3f   "%state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    for k, v in getEnergyDecomp(simulation.context, fgrps).items():                            
        forcenames+="%s "%k.__class__.__name__[:-5]
        values+=str("%.3f   "%v.value_in_unit(kilojoules_per_mole))
    if printHeaders:
        print(forcenames)
    print(values)


def getEnergyDecomp(context, forcegroups):
    energies = {}
    for f, i in forcegroups.items():
        energies[f] = context.getState(getEnergy=True,
                                       groups=1 << i).getPotentialEnergy()
    return energies

###to compute the off term
def offfunc(d1,d2):
    A1=358.679211
    A2=11860.007062
    r1=0.114303
    r2=0.111936
    D1=1346.049696
    D2=94.610940
    return (D1*np.exp(-A1*(d1-r1)**2)+D2*np.exp(-A2*(d2-r2)**2)) # in kJ/mol

###to compute distance with PBC
@jit(float64(float64[:],float64[:],float64[:]),nopython=True)
def calcPBC_sqdistance(COM1,COM2,box):
    diff = COM1 - COM2
    for d in range(3):
        while diff[d] >= box[d]/2.:
            diff[d] -= box[d]
        while diff[d] < -box[d]/2.:
            diff[d] += box[d]
    return float64(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2])

###to compute distance, no PBC
@jit(float64(float64[:],float64[:]),nopython=True)
def calc_sqdistance(COM1,COM2):
    diff = COM1 - COM2
    return diff[0]**2 + diff[1]**2 + diff[2]**2

###hand calc EVB energy and c
def check_EVB(simulation,simulation1,simulation2,V,box,isPBC,idx1,idx2,idx3):
    state = simulation.context.getState(getEnergy=True,getPositions=True)
    evbpos = state.getPositions(asNumpy=True)

    ###set the position in state1 and state2 simulations
    #simulation1.context.setPositions(evbpos)
    #simulation2.context.setPositions(evbpos)
    state1 = simulation1.context.getState(getEnergy=True,getPositions=True)
    state2 = simulation2.context.getState(getEnergy=True,getPositions=True)

    EEVB = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole) ###EVB energy
    E1 = state1.getPotentialEnergy().value_in_unit(kilojoules_per_mole) ###state1 PE
    E2 = state2.getPotentialEnergy().value_in_unit(kilojoules_per_mole)

    #print (EEVB, E1, E2)

    #pos1 = state1.getPositions(asNumpy=True)
    #pos2 = state2.getPositions(asNumpy=True)
    #print (pos1[0:5])
    #print (pos2[0:5])

    ###get E12
    pos12 = state.getPositions(asNumpy=True).value_in_unit(nanometer)
    if isPBC:
        d1 = np.sqrt(calcPBC_sqdistance(pos12[idx1],pos12[idx2],box)) ##get distance between O--H
        d2 = np.sqrt(calcPBC_sqdistance(pos12[idx2],pos12[idx3],box)) ##get distance between H--C
    else:
        d1 = np.sqrt(calc_sqdistance(pos12[idx1],pos12[idx2])) ##get distance between O--H
        d2 = np.sqrt(calc_sqdistance(pos12[idx2],pos12[idx3])) ##get distance between H--C
    E12 = offfunc(d1,d2) ##compute the off term by hand

    ###
    EEVB = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole) ###EVB energy
    E1 = state1.getPotentialEnergy().value_in_unit(kilojoules_per_mole) ###state1 PE
    E2 = state2.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    E2V = E2 + V ###state2 PE
    EVB1 = 0.5*(E1 + E2V - np.sqrt((E1-E2V)**2+4*E12**2))

    T = E1 + E2V
    D = E1*E2V - E12*E12
    eigenvalue1 = (T + np.sqrt((E1-E2V)**2+4*E12**2)) * 0.5
    eigenvalue2 = (T - np.sqrt((E1-E2V)**2+4*E12**2)) * 0.5

    #print ("##eigen: ", eigenvalue1, eigenvalue2)

    mineigen = min(eigenvalue1,eigenvalue2)

    #print ("##mineigen: ", mineigen)

    c = [-E12/(E1-mineigen), 1]
    L = np.sqrt(c[0]**2 + c[1]**2)
    c /= L

    EVB2 = c[0]**2*E1+c[1]**2*E2V+2*c[0]*c[1]*E12 ##myEVB1 is always equal to myEVB2

    
    print ("OMEVB: ", "%f"%(EEVB))
    print ("myEVB1: ", "%f"%(EVB1))
    print ("myEVB2: ", "%f"%(EVB2))
    print ("E1: ", "%f"%(E1))
    print ("E2: ", "%f"%(E2))
    print ("E2V: ", "%f"%(E2V))
    print ("E12: ", "%f"%(E12))
    print ("c1^2: %f, c1: %f"%(c[0]**2,c[0]))
    print ("c2^2: %f, c2: %f"%(c[1]**2,c[1]))
    print ("d1: ", d1)                                                                                            
    print ("d2: ", d2)  

    return EEVB, EVB1, EVB2, E1, E2, E12, c[0], c[1], d1, d2

def newswap_idx(simulation,swappedSimulation,box,isPBC,Cidx,watHidx,iswater):
    
    swapped = False
    


    state = simulation.context.getState(getPositions=True,getVelocities=True,getEnergy=True)    
    pos = state.getPositions(asNumpy=True)._value    
    vec = state.getVelocities(asNumpy=True)


    if isPBC: distsq = calcPBC_sqdistance(pos[0],pos[4],box)
    else: distsq = calc_sqdistance(pos[0],pos[4])
    #print("dist between 0 and 4 %f"%np.sqrt(distsq))




    if iswater: #c2 > c1
        #print ("it is water")
        Cpos = pos[Cidx]
        Hpos = pos[watHidx]

        distsq = np.zeros(len(Hpos))
        for i in range(len(Hpos)):
            if isPBC: distsq[i] = calcPBC_sqdistance(Cpos,Hpos[i],box)
            else: distsq[i] = calc_sqdistance(Cpos,Hpos[i])

        idx = int(np.argmin(distsq)/2) ## get the id of the closest water
        hidx = np.argmin(distsq)%2 ## determine which hydrogen in water is the closest

        if idx != 0: ##swap water
            startidx = (idx - 1) * 5 + 11 ##start from 2nd water O
            endidx = (idx - 1) * 5 + 16 ##ends at 2nd water H2, + 1 for range()
            pos[[range(0,5),range(startidx,endidx)]] = pos[[range(startidx,endidx),range(0,5)]]
            vec[[range(0,5),range(startidx,endidx)]] = vec[[range(startidx,endidx),range(0,5)]]
            
         #   print ("swap water proposal")
            swapped = True

        if hidx%2 == 0: ##swap hydrogen position if the 1st is the closest
            pos[[3,4]] = pos[[4,3]]
            vec[[3,4]] = vec[[4,3]]
            swapped = True
          #  print ("swap hydrogen proposal")
        wrappos(state,box,range(5),5)

    ### if it is OHR:
    else:
        print ("it is OHR")
        CHpos = pos[[7,8]]
        Opos = pos[0]
        distsq = [distsq,0,0]
        if isPBC: 
            distsq[1] = calcPBC_sqdistance(Opos,CHpos[0],box)
            distsq[2] = calcPBC_sqdistance(Opos,CHpos[1],box)
        else: 
            distsq[1] = calc_sqdistance(Opos,CHpos[0])
            distsq[2] = calc_sqdistance(Opos,CHpos[1])
        idx = int(np.argmin(distsq))
        if idx != 0: ##closest H in CH3 from oxygen of first water has changed

            pos[[4,6+idx]] = pos[[6+idx,4]]
            vec[[4,6+idx]] = vec[[6+idx,4]]
            swapped = True
           # print ("closest CH3 H change proposal")

    # if swapped:
    #     swappedSimulation.context.setPositions(pos)
    #     swappedSimulation.context.setVelocities(vec)
    #     swappedState = swappedSimulation.context.getState(getEnergy= True)
    #     swapEEVB = swappedState.getPotentialEnergy()._value
    #     EEVB = state.getPotentialEnergy()._value
    #     if swapEEVB < EEVB:
    #         simulation.context.setPositions(pos)
    #         simulation.context.setVelocities(vec)
    #         print("\n\n\nSwap happened\n\n\n")


def swap_idx(simulation,box,isPBC,Cidx,watHidx,iswater):
    state = simulation.context.getState(getPositions=True,getVelocities=True,getEnergy=True)
    
    pos = state.getPositions(asNumpy=True)._value    
    vec = state.getVelocities(asNumpy=True)
        

 

    if isPBC: distsq = calcPBC_sqdistance(pos[0],pos[4],box)
    else: distsq = calc_sqdistance(pos[0],pos[4])
    print("dist between 0 and 4 %f"%np.sqrt(distsq))

    #if distsq < 0.15*0.15: iswater = True
    #else: iswater = False

    ### if it is water
    if iswater: #c2 > c1
        print ("it is water")
        Cpos = pos[Cidx]
        Hpos = pos[watHidx]

        distsq = np.zeros(len(Hpos))
        for i in range(len(Hpos)):
            if isPBC: distsq[i] = calcPBC_sqdistance(Cpos,Hpos[i],box)
            else: distsq[i] = calc_sqdistance(Cpos,Hpos[i])

        idx = int(np.argmin(distsq)/2) ## get the id of the closest water
        hidx = np.argmin(distsq)%2 ## determine which hydrogen in water is the closest

        if idx != 0: ##swap water
            startidx = (idx - 1) * 5 + 11 ##start from 2nd water O
            endidx = (idx - 1) * 5 + 16 ##ends at 2nd water H2, + 1 for range()
            pos[[range(0,5),range(startidx,endidx)]] = pos[[range(startidx,endidx),range(0,5)]]
            vec[[range(0,5),range(startidx,endidx)]] = vec[[range(startidx,endidx),range(0,5)]]
            
            print ("swap water")

        if hidx%2 == 0: ##swap hydrogen position if the 1st is the closest
            pos[[3,4]] = pos[[4,3]]
            vec[[3,4]] = vec[[4,3]]
            print ("swap hydrogen")
        wrappos(state,box,range(5),5)

    ### if it is OHR:
    else:
        print ("it is OHR")
        CHpos = pos[[7,8]]
        Opos = pos[0]
        distsq = [distsq,0,0]
        if isPBC: 
            distsq[1] = calcPBC_sqdistance(Opos,CHpos[0],box)
            distsq[2] = calcPBC_sqdistance(Opos,CHpos[1],box)
        else: 
            distsq[1] = calc_sqdistance(Opos,CHpos[0])
            distsq[2] = calc_sqdistance(Opos,CHpos[1])
        idx = int(np.argmin(distsq))
        if idx != 0: ##closest H in CH3 from oxygen of first water has changed

            pos[[4,6+idx]] = pos[[6+idx,4]]
            vec[[4,6+idx]] = vec[[6+idx,4]]
            print ("closest CH3 H changed")

    simulation.context.setPositions(pos)
    simulation.context.setVelocities(vec)
 
