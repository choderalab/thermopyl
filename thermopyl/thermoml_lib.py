import re
import copy
import pandas as pd
import thermoml_schema  # Obtained by `wget http://media.iupac.org/namespaces/ThermoML/ThermoML.xsd` and `pyxbgen ThermoML.xsd`

class Parser(object):
    def __init__(self, filename):
        """Create a parser object from an XML filename."""
        self.filename = filename
        self.root = thermoml_schema.CreateFromDocument(open(self.filename).read())
        
        self.store_compounds()
    
    def store_compounds(self):
        """Extract and store compounds from a thermoml XML file."""
        self.compound_num_to_name = {}
        self.compound_name_to_formula = {}
        for Compound in self.root.Compound:
            nOrgNum = Compound.RegNum.nOrgNum
            sCommonName = Compound.sCommonName[0]
            sFormulaMolec = Compound.sFormulaMolec

            self.compound_num_to_name[nOrgNum] = sCommonName
            self.compound_name_to_formula[sCommonName] = sFormulaMolec

    def parse(self):
        """Parse the current XML filename and return a list of measurements."""
        alldata = []
        for PureOrMixtureData in self.root.PureOrMixtureData:
            components = []
            for Component in PureOrMixtureData.Component:
                nSampleNm = Component.nSampleNm
                nOrgNum = Component.RegNum.nOrgNum
                sCommonName = self.compound_num_to_name[nOrgNum]
                components.append(sCommonName)

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


def count_atoms(formula_string):
    """Parse a chemical formula and return the total number of atoms."""
    element_counts = formula_to_element_counts(formula_string)
    return sum(val for key, val in element_counts.items())


def count_atoms_in_set(formula_string, which_atoms):
    """Parse a chemical formula and return the number of atoms in a set of atoms."""
    element_counts = formula_to_element_counts(formula_string)
    return sum(val for key, val in element_counts.items() if key in which_atoms)


def get_first_entry(cas):
    """If cirpy returns several CAS results, extracts the first one."""
    if type(cas) == type([]):
        cas = cas[0]
    return cas


def formula_to_element_counts(formula_string):
    """Transform a chemical formula into a dictionary of (element, number) pairs."""
    pattern = r'([A-Z][a-z]{0,2}\d*)'
    pieces = re.split(pattern, formula_string)
    #print "\nformula_string=%r pieces=%r" % (formula_string, pieces)
    data = pieces[1::2]
    rubbish = filter(None, pieces[0::2])
    pattern2 = r'([A-Z][a-z]{0,2})'

    results = {}
    for piece in data:
        #print(piece)
        element, number = re.split(pattern2, piece)[1:]
        try:
            number = int(number)
        except ValueError:
            number = 1
        results[element] = number

    return results
