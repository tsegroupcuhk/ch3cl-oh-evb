import copy
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import evb.utils as utils

class EVBSystem:
    def __init__(self, sys1, sys2, sysoff, isPBC):
        self.sub1 = sys1
        self.sub2 = sys2
        self.suboff = sysoff
        self.boxv = sys1.getDefaultPeriodicBoxVectors()
        self.isPBC = isPBC

    def mergeSystem(self, V):
        sysall = mm.System()
        for n in range(self.sub2.getNumParticles()):
            sysall.addParticle(self.sub2.getParticleMass(n))
        # set virtual sites
        for nparticle in range(self.sub2.getNumParticles()):
            if self.sub2.isVirtualSite(nparticle):
                vs = self.sub2.getVirtualSite(nparticle)
                vsnp = vs.getNumParticles()
                vsw = [vs.getWeight(i) for i in range(vsnp)]
                vsid = [vs.getParticle(i) for i in range(vsnp)]
                if vsnp == 2:
                    newvs = mm.TwoParticleAverageSite(vsid[0],vsid[1],vsw[0],vsw[1])
                elif vsnp == 3:
                    newvs = mm.ThreeParticleAverageSite(vsid[0],vsid[1],vsid[2],vsw[0],vsw[1],vsw[2])
                sysall.setVirtualSite(nparticle, newvs)
        # set constraints
        for ncons in range(self.sub2.getNumConstraints()):
            p1, p2, dist = self.sub2.getConstraintParameters(ncons)
            sysall.addConstraint(p1, p2, dist)

        newforces1 = {}
        newdrude = None
        forcecount1 = {}
        for nf, force in enumerate(self.sub1.getForces()):
            forcename = force.__class__.__name__
            forcecount1[forcename] = 0
        for nf, force in enumerate(self.sub1.getForces()):
            forcename = force.__class__.__name__
            forcecount1[forcename] += 1
            newforces1["%s%d"%(forcename,forcecount1[forcename])] = copy.deepcopy(force)
            if isinstance(force, mm.DrudeForce):
                newdrude = mm.DrudeForce()
                for ndrude in range(force.getNumParticles()):
                    p, p1, p2, p3, p4, chrg, polar, a12, a34 = force.getParticleParameters(ndrude)
                    newdrude.addParticle(p, p1, p2, p3, p4, 0.0, polar, a12,a34)

        energy = "+".join([_ for _ in newforces1.keys()])
        sys1f = mm.CustomCVForce(energy)
        for k, v in newforces1.items():
            sys1f.addCollectiveVariable(k, v)

        #print (energy)
        #exit()

        newforces2 = {}
        forcecount2 = {}
        for nf, force in enumerate(self.sub2.getForces()):
            forcename = force.__class__.__name__
            forcecount2[forcename] = 0
        for nf, force in enumerate(self.sub2.getForces()):
            forcename = force.__class__.__name__
            forcecount2[forcename] += 1
            newforces2["%s%d"%(forcename,forcecount2[forcename])] = copy.deepcopy(force)
        energy = "+".join([_ for _ in newforces2.keys()])
        sys2f = mm.CustomCVForce(energy)
        for k, v in newforces2.items():
            sys2f.addCollectiveVariable(k, v)

        newforcesoff = {}
        for nf, force in enumerate(self.suboff.getForces()):
            #print(force)
            if isinstance(force, mm.CustomCompoundBondForce):
                newforcesoff["eneroff"] = copy.deepcopy(force)
            else:
                #print(force, "N")
                pass
        energy = "+".join([_ for _ in newforcesoff.keys()])
        foff = mm.CustomCVForce(energy)
        for k, v in newforcesoff.items():
            foff.addCollectiveVariable(k, v)

        sysf = mm.CustomCVForce(
            "0.5*(E1+E2v-sqrt(DE^2+4*E12^2)); DE=E1-E2v; E2v=E2+V")
        sysf.addCollectiveVariable("E1", sys1f)
        sysf.addCollectiveVariable("E2", sys2f)
        sysf.addCollectiveVariable("E12", foff)
        sysf.addGlobalParameter("V", V)
        sysall.addForce(sysf)

        if newdrude is not None:
            sysall.addForce(newdrude)

        comremover = mm.CMMotionRemover()
        comremover.setFrequency(1000)
        sysall.addForce(comremover)
        if self.isPBC:
            sysall.setDefaultPeriodicBoxVectors(self.boxv[0],self.boxv[1],self.boxv[2])
        return sysall
