import openmm as mm
import openmm.app as app
import numpy as np
from openmm.unit import *
from numba import *
from evb.utils import *
###copy this file to the evb folder

### for each pair of atomtypes A, B
def nbfix_type_decomposition(system,simulation):
    nbfix = [system.getForce(i) for i in range(system.getNumForces()) if isinstance(system.getForce(i), mm.CustomNonbondedForce)][0]

    numAtomClasses = nbfix.getTabulatedFunction(0).getFunctionParameters()[0]
    epstable = nbfix.getTabulatedFunction(0).getFunctionParameters()[2]
    sigtable = nbfix.getTabulatedFunction(1).getFunctionParameters()[2]

    print ("##customnonbonded type decomposition##")

    mysumPE = 0
    for i in range(numAtomClasses):
        for j in range(i,numAtomClasses):
            ## get the target pair parameter
            paireps = epstable[i*numAtomClasses+j]
            ## create a new table
            newtable = np.zeros(numAtomClasses*numAtomClasses)
            ## set the target pair parameter
            newtable[i*numAtomClasses+j] = paireps
            newtable[j*numAtomClasses+i] = paireps
            ## update the table within nbfix
            nbfix.getTabulatedFunction(0).setFunctionParameters(numAtomClasses,numAtomClasses,newtable)
            ## update nbfix within simulation.context
            nbfix.updateParametersInContext(simulation.context)

            ### get the energies
            mystate = simulation.context.getState(getEnergy=True)
            myfgrps = forcegroupify(system)
            for k, v in getEnergyDecomp(simulation.context, myfgrps).items():
                if k.__class__.__name__ == "CustomNonbondedForce":
                    PE = v.value_in_unit(kilojoules_per_mole)
                    print(i, j, PE)
                    mysumPE += PE
    print ("nbfix: ", mysumPE)

    ## reset to the original table
    nbfix.getTabulatedFunction(0).setFunctionParameters(numAtomClasses,numAtomClasses,epstable)
    nbfix.updateParametersInContext(simulation.context)

### for each atom pair indexes i, j
def nbfix_pair_decomposition(system,simulation):
    nbfix = [system.getForce(i) for i in range(system.getNumForces()) if isinstance(system.getForce(i), mm.CustomNonbondedForce)][0]

    numAtomClasses = nbfix.getTabulatedFunction(0).getFunctionParameters()[0]
    epstable = nbfix.getTabulatedFunction(0).getFunctionParameters()[2]
    sigtable = nbfix.getTabulatedFunction(1).getFunctionParameters()[2]

    print ("##customnonbonded pair decomposition##")

    mysumPE = 0
    for i in range(numAtomClasses):
        paireps = epstable[i*numAtomClasses+i]
        if paireps == 0.0:
           zeroPEidx = i
           break

    oridx = np.zeros(nbfix.getNumParticles())
    ##get original atom types:
    for i in range(nbfix.getNumParticles()):
        oridx[i] = nbfix.getParticleParameters(i)[0]
        ## turn all pair to zeroPEidx to turn them off
        nbfix.setParticleParameters(i,[zeroPEidx])

    for i in range(nbfix.getNumParticles()):
        for j in range(i,nbfix.getNumParticles()):
            ##only turn on one pair interaction
            nbfix.setParticleParameters(i,[oridx[i]])
            nbfix.setParticleParameters(j,[oridx[j]])

            ## update nbfix within simulation.context
            nbfix.updateParametersInContext(simulation.context)

            ### get the energies
            mystate = simulation.context.getState(getEnergy=True)
            myfgrps = forcegroupify(system)
            for k, v in getEnergyDecomp(simulation.context, myfgrps).items():
                if k.__class__.__name__ == "CustomNonbondedForce":
                    PE = v.value_in_unit(kilojoules_per_mole)
                    print(i, j, PE)
                    mysumPE += PE

            ## turn off again
            nbfix.setParticleParameters(i,[zeroPEidx])
            nbfix.setParticleParameters(j,[zeroPEidx])

    print ("nbfix: ", mysumPE)

    ## reset to the original idx
    for i in range(nbfix.getNumParticles()):
        nbfix.setParticleParameters(i,[oridx[i]])
    nbfix.updateParametersInContext(simulation.context)

###does not work for PBC systems
def get_coul_pair_decomposition(system,simulation):
    for i in range(system.getNumForces()):
        sysforce = system.getForce(i)
        if isinstance(sysforce, mm.NonbondedForce):
            #print ("Use Dispersion Correction: ", sysforce.getUseDispersionCorrection())
            print ("No. of particles: ", sysforce.getNumParticles())
            print ("No. of exceptions: ", sysforce.getNumExceptions())

            ### check exceptions 1-4 parameters
            firstnonzero = True
            for exidx in range(sysforce.getNumExceptions()):
                p1, p2, qq, sig, eps = sysforce.getExceptionParameters(exidx)

                ### exception parameters cannot be changed 
                ### without creating a new simulation.context
                ### so the 1-4 energy is subtracted instead
                if qq.value_in_unit(elementary_charge**2) != 0.0 or eps.value_in_unit(kilojoules_per_mole) != 0.0:
                    if firstnonzero:
                        print ("there are non-zero 1-4 parameters (will be subtracted from total)")
                        firstnonzero = False
                    print (p1,p2,qq,sig,eps)

            ### store the original parameters
            chargeL = []

            ### get the original parameters
            for atom1 in range(sysforce.getNumParticles()):
                charge, _, _ = sysforce.getParticleParameters(atom1)
                chargeL.append(charge)
                ### set the parameters to zero eps and charge, 1.0 sigma
                sysforce.setParticleParameters(atom1,0.0,1.0,0.0)

            ##check if PE from 1-4
            sysforce.updateParametersInContext(simulation.context)
            mystate = simulation.context.getState(getEnergy=True)
            myfgrps = forcegroupify(system)
            for k, v in getEnergyDecomp(simulation.context, myfgrps).items():
                if k.__class__.__name__ == "NonbondedForce":
                    PE14 = v.value_in_unit(kilojoules_per_mole)

            print ("###Coul decomposition###")

            mysumPE = 0
            ### get contribution of each pair
            for atom1 in range(sysforce.getNumParticles()):
                for atom2 in range(atom1,sysforce.getNumParticles()):
                    sysforce.setParticleParameters(atom1,chargeL[atom1],1.0,0.0)
                    sysforce.setParticleParameters(atom2,chargeL[atom2],1.0,0.0)
                    
                    ## update sysforce within simulation.context
                    sysforce.updateParametersInContext(simulation.context)

                    ### get the energies
                    mystate = simulation.context.getState(getEnergy=True)
                    myfgrps = forcegroupify(system)
                    for k, v in getEnergyDecomp(simulation.context, myfgrps).items():
                        if k.__class__.__name__ == "NonbondedForce":
                            PE = v.value_in_unit(kilojoules_per_mole) - PE14
                            print(atom1, atom2, PE)
                            mysumPE += PE
                    ## reset charge to 0
                    sysforce.setParticleParameters(atom1,0.0,1.0,0.0)
                    sysforce.setParticleParameters(atom2,0.0,1.0,0.0)
            print ("nonbonded(Coul): ", mysumPE)
            print ("1-4PE: ", PE14)
            print ("nonbonded(total): ", mysumPE + PE14)

            ## reset parameters:
            for atom1 in range(sysforce.getNumParticles()):
                sysforce.setParticleParameters(atom1,chargeL[atom1],1.0,0.0)
            sysforce.updateParametersInContext(simulation.context)
            break

def check_harmonic_bond_energy(simulation,system,box,isPBC):
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
    sumPE = 0
    for i in range(system.getNumForces()):
        sysforce = system.getForce(i)
        if isinstance(sysforce, mm.HarmonicBondForce):
            for j in range(sysforce.getNumBonds()):
                ll = sysforce.getBondParameters(j)
                index1 = ll[0]
                index2 = ll[1]
                r0 = ll[2].value_in_unit(nanometer)
                if isPBC: r = np.sqrt(calcPBC_sqdistance(pos[index1],pos[index2],box))
                else: r = np.sqrt(calc_sqdistance(pos[index1],pos[index2]))
                k = sysforce.getBondParameters(j)[3].value_in_unit(kilojoule/(nanometer**2*mole))
                bondPE = 0.5*k*(r-r0)**2
                #print("bond length %.4f bond energy %.4f"%(r,bondPE))
                sumPE += bondPE
            print("harmonicbondPE: ", sumPE)
            break

def check_drude_energy(simulation,system,box,isPBC):
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
    sumPE = 0
    energydipole = 0
    ONE_4PI_EPS0 = 138.935456 ## kilojoules_per_mole*nanometer

    for i in range(system.getNumForces()):
        sysforce = system.getForce(i)
        if isinstance(sysforce, mm.DrudeForce):
            for j in range(sysforce.getNumParticles()):
                ll = sysforce.getParticleParameters(j)
                index1=ll[0]
                index2=ll[1]
                charge=ll[5]
                if isPBC: r = np.sqrt(calcPBC_sqdistance(pos[index1],pos[index2],box))
                else: r = np.sqrt(calc_sqdistance(pos[index1],pos[index2]))
                k = 418400.0 #kJ/mol nm^2
                bondPE = 0.5*k*r**2
                print("\natom pair %d %d bond length %.4f bond energy %.4f"%(index1,index2,r,bondPE),'charge',charge)
                print(' coord:',pos[index1],pos[index2]) 
                sumPE += bondPE
         
            print("\nDrudePE: ", sumPE)
            break

def computeScreening(r, sumthole, alpha1, alpha2):
    u = r * sumthole / (alpha1*alpha2)**(1.0/6.0)
    return 1.0 - (1.0+u/2)*np.exp(-u)

def check_drude_energy_screen(simulation,system,box,isPBC):
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
    sumPE = 0
    energydipole = 0
    ONE_4PI_EPS0 = 138.935456 ## kilojoules_per_mole*nanometer

    ##get nonbonded force for charges
    for i in range(system.getNumForces()):
        nbforce = system.getForce(i)
        if isinstance(nbforce, mm.NonbondedForce):
            break

    for i in range(system.getNumForces()):
        sysforce = system.getForce(i)
        if isinstance(sysforce, mm.DrudeForce):
            for j in range(sysforce.getNumParticles()):
                ll = sysforce.getParticleParameters(j)
                index1=ll[0]
                index2=ll[1]
                if isPBC: r = np.sqrt(calcPBC_sqdistance(pos[index1],pos[index2],box))
                else: r = np.sqrt(calc_sqdistance(pos[index1],pos[index2]))
                k = 418400.0 #kJ/mol nm^2
                bondPE = 0.5*k*r**2
                #print("bond length %.4f bond energy %.4f"%(r,bondPE))
                sumPE += bondPE
            for j in range(sysforce.getNumScreenedPairs()):
                p1, p2, sumthole = sysforce.getScreenedPairParameters(j)
                #sysforce.setScreenedPairParameters(j,p1,p2,0.0)
                #sysforce.updateParametersInContext(simulation.context)

                ##the index refers to the index within the defined force
                ##not the system atom indexes
                ##here the index means the x-th atom with drude
                print (p1, p2, sumthole)

                a1, a2, _, _, _, _, alpha1, _, _ = sysforce.getParticleParameters(p1)
                b1, b2, _, _, _, _, alpha2, _, _ = sysforce.getParticleParameters(p2)

                alpha1 = alpha1.value_in_unit(nanometer**3)
                alpha2 = alpha2.value_in_unit(nanometer**3)

                qa1, _, _ = nbforce.getParticleParameters(a1)
                qa2, _, _ = nbforce.getParticleParameters(a2)
                qb1, _, _ = nbforce.getParticleParameters(b1)
                qb2, _, _ = nbforce.getParticleParameters(b2)

                qa1 = qa1.value_in_unit(elementary_charge)
                qa2 = qa2.value_in_unit(elementary_charge)
                qb1 = qb1.value_in_unit(elementary_charge)
                qb2 = qb2.value_in_unit(elementary_charge)

                print (qa1, qa2, qb1, qb2)

                if isPBC:
                    r11 = np.sqrt(calcPBC_sqdistance(pos[a1],pos[b1],box))
                    r12 = np.sqrt(calcPBC_sqdistance(pos[a1],pos[b2],box))
                    r21 = np.sqrt(calcPBC_sqdistance(pos[a2],pos[b1],box))
                    r22 = np.sqrt(calcPBC_sqdistance(pos[a2],pos[b2],box))
                else:
                    r11 = np.sqrt(calc_sqdistance(pos[a1],pos[b1]))
                    r12 = np.sqrt(calc_sqdistance(pos[a1],pos[b2]))
                    r21 = np.sqrt(calc_sqdistance(pos[a2],pos[b1]))
                    r22 = np.sqrt(calc_sqdistance(pos[a2],pos[b2]))

                energydipole += ONE_4PI_EPS0*qa1*qb1*computeScreening(r11,sumthole,alpha1,alpha2)/r11
                energydipole += ONE_4PI_EPS0*qa1*qb2*computeScreening(r12,sumthole,alpha1,alpha2)/r12
                energydipole += ONE_4PI_EPS0*qa2*qb1*computeScreening(r21,sumthole,alpha1,alpha2)/r21
                energydipole += ONE_4PI_EPS0*qa2*qb2*computeScreening(r22,sumthole,alpha1,alpha2)/r22

            print("sumDrudePE: ", sumPE, " energydipole: ", energydipole, " added: ", sumPE+energydipole)
            break


