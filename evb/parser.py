import copy
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np


class SN2OffDiagGenerator:
    def __init__(self, hamiltonian):
        self.ff = hamiltonian
        self.params = {}
        self._jaxPotential = None
        self.types = []

    def registerPairType(self, angle):
        types = self.ff._findAtomTypes(angle, 3)
        self.types = types
        self.params['D'] = float(angle['D'])
        self.params['A'] = float(angle['A'])
        self.params['d0'] = float(angle['d0'])
        self.params['a'] = float(angle['a'])
        self.params['b'] = float(angle['b'])
        self.params['c'] = float(angle['c'])
        

    @staticmethod
    def parseElement(element, hamiltonian):
        generator = SN2OffDiagGenerator(hamiltonian)
        hamiltonian.registerGenerator(generator)
        for bondtype in element.findall("Pair"):
            generator.registerPairType(bondtype.attrib)

    def createForce(self, system, data, nonbondedMethod, nonbondedCutoff,
                    args):

        n_angles = len(data.angles)
        # build map
        ifFound = False
        for i in range(n_angles):
            idx1 = data.angles[i][0]
            idx2 = data.angles[i][1]
            idx3 = data.angles[i][2]
            type1 = data.atomType[data.atoms[idx1]]
            type2 = data.atomType[data.atoms[idx2]]
            type3 = data.atomType[data.atoms[idx3]]
            if type2 in self.types[1]:
                if type1 in self.types[0] and type3 in self.types[2]:
                    ifFound = True
                    atom1 = idx1
                    atom2 = idx2
                    atom3 = idx3
                elif type1 in self.types[2] and type3 in self.types[0]:
                    ifFound = True
                    atom1 = idx3
                    atom2 = idx2
                    atom3 = idx1
            if ifFound:
                break
        if not ifFound:
            raise BaseException("No pair found for EVB off diagonal term")

        force = mm.CustomCompoundBondForce(
            3,
            "D*exp(-A*(d2-d1-d0)^2)*(a*(d2-d1)^2+b*(d2-d1)+c); d1=distance(p1,p2); d2=distance(p2,p3)"
        )
        force.addPerBondParameter("D")        
        force.addPerBondParameter("A")
        force.addPerBondParameter("d0")
        force.addPerBondParameter("a")
        force.addPerBondParameter("b")
        force.addPerBondParameter("c")
        force.addBond([atom1, atom2, atom3], [
            self.params["D"],self.params["A"],self.params["d0"], self.params["a"], self.params["b"],
            self.params["c"]
        ])
        system.addForce(force)


app.forcefield.parsers["SN2OffDiagForce"] = SN2OffDiagGenerator.parseElement
