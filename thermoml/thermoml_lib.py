import copy
import pandas as pd
import thermoml_schema  # Obtained by `wget http://media.iupac.org/namespaces/ThermoML/ThermoML.xsd` and `pyxbgen ThermoML.xsd`

class Parser(object):
    def __init__(self, filename):
        print(filename)
        self.filename = filename
        self.root = thermoml_schema.CreateFromDocument(open(self.filename).read())
        
        self.store_compounds()
    
    def store_compounds(self):
        self.compound_num_to_name = {}
        self.compound_name_to_formula = {}
        for Compound in self.root.Compound:
            nOrgNum = Compound.RegNum.nOrgNum
            sCommonName = Compound.sCommonName[0]
            sFormulaMolec = Compound.sFormulaMolec
            
            self.compound_num_to_name[nOrgNum] = sCommonName
            self.compound_name_to_formula[sCommonName] = sFormulaMolec

    def parse(self):        
        alldata = []
        for PureOrMixtureData in self.root.PureOrMixtureData:
            components = []
            for Component in PureOrMixtureData.Component:
                nSampleNm = Component.nSampleNm
                nOrgNum = Component.RegNum.nOrgNum
                sCommonName = self.compound_num_to_name[nOrgNum]
                print(sCommonName)
                components.append(sCommonName)
       
            print(components)
            components_string = "__".join(components)

            property_dict = {}
            for Property in PureOrMixtureData.Property:
                nPropNumber = Property.nPropNumber
                ePropName = Property.Property_MethodID.PropertyGroup.content()[0].ePropName  # ASSUMING LENGTH 1
                property_dict[nPropNumber] = ePropName

            state = dict(filename=self.filename, components=components_string)
            
            state["Pressure, kPa"] = None  # This is the only pressure unit used in ThermoML
            state['Temperature, K'] = None  # This is the only temperature unit used in ThermoML
            
            composition = dict()
            for Constraint in PureOrMixtureData.Constraint:
                nConstraintValue = Constraint.nConstraintValue
                ConstraintType = Constraint.ConstraintID.ConstraintType
                
                assert len(ConstraintType.content()) == 1
                constraint_type = ConstraintType.content()[0]
                state[constraint_type] = nConstraintValue
                
                if constraint_type in ["Mole fraction", "Mass Fraction", "Molality, mol/kg", "Solvent: Amount concentration (molarity), mol/dm3"]:
                    nOrgNum = Constraint.ConstraintID.RegNum.nOrgNum
                    sCommonName = self.compound_num_to_name[nOrgNum]
                    if Constraint.Solvent is not None:
                        solvents = [self.compound_num_to_name[x.nOrgNum] for x in Constraint.Solvent.RegNum]
                    else:
                        solvents = []                
                    solvent_string = "%s___%s" % (sCommonName, "__".join(solvents))
                    state["%s metadata" % constraint_type] = solvent_string

            variable_dict = {}
            for Variable in PureOrMixtureData.Variable:
                nVarNumber = Variable.nVarNumber
                VariableType = Variable.VariableID.VariableType
                assert len(VariableType.content()) == 1
                vtype = VariableType.content()[0]  # Assume length 1, haven't found counterexample yet.
                variable_dict[nVarNumber] = vtype
                if vtype in ["Mole fraction", "Mass Fraction", "Molality, mol/kg", "Solvent: Amount concentration (molarity), mol/dm3"]:
                    nOrgNum = Variable.VariableID.RegNum.nOrgNum
                    sCommonName = self.compound_num_to_name[nOrgNum]
                    if Variable.Solvent is not None:
                        solvents = [self.compound_num_to_name[x.nOrgNum] for x in Variable.Solvent.RegNum]
                    else:
                        solvents = []
                    solvent_string = "%s___%s" % (sCommonName, "__".join(solvents))
                    state["%s Variable metadata" % vtype] = solvent_string
            
            
            for NumValues in PureOrMixtureData.NumValues:
                current_data = copy.deepcopy(state)  # Copy in values of constraints.
                current_composition = copy.deepcopy(composition)
                for VariableValue in NumValues.VariableValue:
                    nVarValue = VariableValue.nVarValue
                    nVarNumber = VariableValue.nVarNumber
                    vtype = variable_dict[nVarNumber]
                    current_data[vtype] = nVarValue

                for PropertyValue in NumValues.PropertyValue:
                    nPropNumber = PropertyValue.nPropNumber
                    nPropValue = PropertyValue.nPropValue
                    ptype = property_dict[nPropNumber]
                    current_data[ptype] = nPropValue

                alldata.append(current_data)
        return alldata


def is_good(formula_string):
    good = set(["H", "N", "C", "O", "S", ""])
    elements = set(re.findall("[A-Z][a-z]?", formula_string))    
    if good.intersection(elements) is None:  # Has no good elements
        return False
    if len(elements.difference(good)) > 0:  # Has an unwanted element
        return False
    return True

def count_heavy_atoms(formula_string):
    heavy_atoms = ["N", "C", "O", "S"]
    elements = re.findall("[A-Z][a-z]?\d?\d?\d?", formula_string)
    print(elements)
    n_heavy = 0
    for s in elements:
        try:
            n_atoms = int(re.split("[A-Z][a-z]?", s)[1])
        except:
            n_atoms = 1
        atom = re.split("\d?\d?\d?", s)[0]
        print(s, atom, n_atoms)
        if atom in heavy_atoms:
            n_heavy += n_atoms
    return n_heavy
        

def get_first_entry(cas):
    """If cirpy returns several CAS results, extracts the first one."""
    if type(cas) == type([]):
        cas = cas[0]
    return cas
