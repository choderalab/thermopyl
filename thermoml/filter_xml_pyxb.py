import copy
import pandas as pd
import glob
import thermoml_schema  # Obtained by `wget http://media.iupac.org/namespaces/ThermoML/ThermoML.xsd` and `pyxbgen ThermoML.xsd`

def parse(filename):
    print(filename)
    root = thermoml_schema.CreateFromDocument(open(filename).read())

    compound_dict = {}
    for Compound in root.Compound:
        nOrgNum = Compound.RegNum.nOrgNum
        sCommonName = Compound.sCommonName[0]
        sFormulaMolec = Compound.sFormulaMolec
        compound_dict[nOrgNum] = dict(sCommonName=sCommonName, sFormulaMolec=sFormulaMolec)

    alldata = []
    for PureOrMixtureData in root.PureOrMixtureData:
        components = []
        for Component in PureOrMixtureData.Component:
            nSampleNm = Component.nSampleNm
            nOrgNum = Component.RegNum.nOrgNum
            sCommonName = compound_dict[nOrgNum]["sCommonName"]
            print(sCommonName)
            #component_list.append(dict(sCommonName=sCommonName, nOrgNum=nOrgNum))            
            components.append(sCommonName)
   
        print(components)
        components_string = "__".join(components)

        property_dict = {}
        for Property in PureOrMixtureData.Property:
            nPropNumber = Property.nPropNumber
            ePropName = Property.Property_MethodID.PropertyGroup.content()[0].ePropName  # ASSUMING LENGTH 1
            property_dict[nPropNumber] = ePropName

        state = dict(filename=filename, components=components_string)
        
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
                sCommonName = compound_dict[nOrgNum]["sCommonName"]
                if Constraint.Solvent is not None:
                    solvents = [compound_dict[x.nOrgNum]["sCommonName"] for x in Constraint.Solvent.RegNum]
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
                sCommonName = compound_dict[nOrgNum]["sCommonName"]
                if Variable.Solvent is not None:
                    solvents = [compound_dict[x.nOrgNum]["sCommonName"] for x in Variable.Solvent.RegNum]
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
    return alldata, root

#filename = "./10.1007/s10765-010-0742-8.xml"
#filename = "./10.1021/je300826t.xml"
#alldata, root = parse(filename)

data = []
for filename in glob.glob("./*/*.xml")[0:250]:
    try:
        alldata, root = parse(filename)
    except IOError:
        continue
    for d in alldata:
        if True or u'Mass density, kg/m3' in d:
            data.append(d)

data = pd.DataFrame(data)
#data.to_hdf("./data.h5", 'data')

X = data.ix[data["Mass density, kg/m3"].dropna().index]
X = X[X["Temperature, K"] > 270]
X = X[X["Temperature, K"] < 330]
X = X[X["Pressure, kPa"] > 50.]
X = X[X["Pressure, kPa"] < 150.]
X.dropna(axis=1, how='all', inplace=True)


chemicals = X.components.apply(lambda x: x.split("__")).values
s = set()
for chemlist in chemicals:
    for chem in chemlist:
        s.add(chem)
chemicals = s
