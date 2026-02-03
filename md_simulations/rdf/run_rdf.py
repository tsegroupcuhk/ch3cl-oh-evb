import openmm as mm
import openmm.app as app
from openmm.unit import *
from sys import stdout, exit, stderr, argv
from datetime import datetime as dtt
import numpy as np
from openmmplumed import *

starttime = dtt.now()
print(starttime)
stdout.flush()

def forcegroupify(system):
    forcegroups = {}
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        force.setForceGroup(i)
        forcegroups[force] = i
    return forcegroups

def getEnergyDecomp(context, forcegroups):
    energies = {}
    for f, i in forcegroups.items():
        energies[f] = context.getState(
            getEnergy=True, groups=1 << i).getPotentialEnergy()
    return energies

def modify_nbfix(system, nbcut, switchdist):
    nbfix = [
        system.getForce(i) for i in range(system.getNumForces())
        if isinstance(system.getForce(i), mm.CustomNonbondedForce)
    ][0]
    nbfix.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    nbfix.setCutoffDistance(nbcut)
    nbfix.setUseSwitchingFunction(True)
    nbfix.setSwitchingDistance(switchdist)

#####################Set steps 
jobname='rdf'
forcefield = app.ForceField('toppar.xml')
Temp=300
nsavcrd=2000 # every 1 ps
nprint=2000 # every 1 ps
nsteps=int(50*2e06) ## 50 ns
dt=0.5*femtoseconds ##use 1.0 fs for DRUDE and flexible H-bonds

# for PBC
ewaldtol = 1e-5
box=np.array([3.104,3.104,3.104])
box_wu = box*nanometer
nbcut = 1.55 * nanometer
switchdist = 1.45 * nanometer

print('###Loading things###')
psf = app.CharmmPsfFile('sys.psf')
crd = app.CharmmCrdFile('sys.crd')

###for nonPBC
# system = forcefield.createSystem(psf.topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
###for PBC
psf.setBox(box_wu[0],box_wu[1],box_wu[2])
system = forcefield.createSystem(psf.topology, nonbondedMethod=app.PME, nonbondedCutoff=nbcut, switchDistance=switchdist, ewaldErrorTolerance=ewaldtol, constraints=None, rigidWater=False)
modify_nbfix(system, nbcut, switchdist)


###plumed script
# RC = distance between O* and O* (O*: oxygen in OH radical)
script = """
UNITS LENGTH=A TIME=fs
dist: DISTANCE ATOMS=5001,5005 NOPBC

restraint: RESTRAINT ARG=dist AT=7 KAPPA=100

PRINT ARG=dist,restraint.bias STRIDE=2000 FILE=rdf.cv
"""
###add plumed force to the system
system.addForce(PlumedForce(script))

###integrator for nonDrude NVT simulations
# integrator = mm.LangevinMiddleIntegrator(Temp*kelvin, 1./picosecond, dt)
# integrator.setConstraintTolerance(1e-06)

### integrator for Drude non-NVE simulations
integrator = mm.DrudeLangevinIntegrator(Temp*kelvin, 5/picosecond, 1*kelvin, 20/picosecond, dt)
integrator.setMaxDrudeDistance(0.02)
integrator.setConstraintTolerance(1e-06)

### integrator for Drude NVE simulations
# integrator = mm.DrudeSCFIntegrator(dt)
# integrator.setMinimizationErrorTolerance(1e-05)

#for energy decomposition later
fgrps = forcegroupify(system)

platname = 'CUDA'
platform = mm.Platform.getPlatformByName(platname)
prop = dict(CudaPrecision='mixed') ###use mixed for NVT, double for NVE
platform.loadPluginsFromDirectory('/home/yutong/software/myOpenMM-plumed-plugin/openmm-plumed1.0_install/lib/plugins/')

simulation = app.Simulation(psf.topology, system, integrator, platform, prop)

###if not read from .rst restart, setup the positions from the .crd and compute virtual sites if any.
simulation.context.setPositions(crd.positions)
simulation.context.computeVirtualSites()
simulation.context.setTime(0)
simulation.context.setStepCount(0)

# print out initial energies
print('#Initial Frame Energies')
box = simulation.context.getState().getPeriodicBoxVectors()
print('Box:', box)
state = simulation.context.getState(getPositions=True, getEnergy=True)
print('PE', state.getPotentialEnergy().value_in_unit(kilojoules_per_mole),'kJ/mol')
for k, v in getEnergyDecomp(simulation.context, fgrps).items():
    print(k.__class__.__name__, v.value_in_unit(kilojoules_per_mole),'kJ/mol')


####set up output files
stdout.flush()
dcd=app.DCDReporter(jobname+'.dcd', nsavcrd)
firstdcdstep = nsavcrd
dcd._dcd = app.DCDFile(dcd._out, simulation.topology,simulation.integrator.getStepSize(), firstdcdstep, nsavcrd)
simulation.reporters.append(dcd)
simulation.reporters.append(app.StateDataReporter(jobname+'.log',nprint, step=True, time=True, kineticEnergy=True, potentialEnergy=True,totalEnergy=True, temperature=True, speed=True,separator=',  '))

###### minimization
#print('###Energy Minimization###')
#stdout.flush()
#simulation.minimizeEnergy() #tolerance=1*kilocalories_per_mole
#state = simulation.context.getState(getPositions=True,getEnergy=True)
#print('PE', state.getPotentialEnergy().value_in_unit(kilojoules_per_mole),'kJ/mol')

## run MD
simulation.step(nsteps)

###create restart file
state = simulation.context.getState(getPositions=True,getVelocities=True,getEnergy=True)
with open(jobname + '.rst', 'w') as f:
    f.write(mm.XmlSerializer.serialize(state))

# print out final energies
print('#Final Frame Energies')
print('PE', state.getPotentialEnergy().value_in_unit(kilojoules_per_mole),'kJ/mol')
for k, v in getEnergyDecomp(simulation.context, fgrps).items():
    print(k.__class__.__name__, v.value_in_unit(kilojoules_per_mole),'kJ/mol')
box = simulation.context.getState().getPeriodicBoxVectors()
print('Box:', box)

endtime = dtt.now()
print(endtime)
print(endtime-starttime)

