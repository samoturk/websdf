from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
import pandas as pd


def _clogSw(mol):
    '''
    Inspired by work by Christos Kannas presented at RDKit UGM 2013
    Based on:
    J. S. Delaney, Journal of Chemical Information and Modeling, 44, 1000-1005,
    2004.
    '''
    MolWeight = Descriptors.MolWt(mol)
    clogP = Descriptors.MolLogP(mol)
    RotBonds = Descriptors.NumRotatableBonds(mol)
    aromaticHeavyatoms = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[a]")))
    numAtoms = mol.GetNumAtoms()
    AromProp = float(aromaticHeavyatoms) / numAtoms
    # New clogSw with coefficients from Christos' presentation
    clogSw_value = 0.233743817233 \
                    -0.74253027 * clogP \
                    -0.00676305 * MolWeight \
                    +0.01580559 * RotBonds \
                    -0.35483585 * AromProp
    return clogSw_value

def _calculate_descs(df, checks):
    # Create a # column
    df['#'] = range(len(df))
    # Put # column on first position
    df = df[['#'] + [x for x in list(df.columns) if x != '#']]    
    # For compatibility with RDKit older than 2015_03
    if 'SMILES' in list(df.columns):
        df = df.drop(['SMILES'], axis=1)
    # Recalculate SMILES
    if 'SMILES' in checks:
        df['SMILES'] = df.apply(lambda x: Chem.MolToSmiles(x['ROMol']), axis=1)
    # Calculate MW
    if 'MW' in checks:
        df['MW'] = df['ROMol'].map(Descriptors.MolWt).round(decimals=2)
    # Calculate logP
    if 'logP' in checks:
        df['logP'] = df['ROMol'].map(Descriptors.MolLogP).round(decimals=2)
    # Calculate H-bond donors and acceptors    
    if 'HB' in checks:
        df['HBA'] = df['ROMol'].map(Descriptors.NumHAcceptors)
        df['HBD'] = df['ROMol'].map(Descriptors.NumHDonors)
    # Calculate solubility    
    if 'logS' in checks:
        df['clogSw'] = df['ROMol'].map(_clogSw).round(decimals=2)
    return df

def read_sdf(sdf, checks):
    df = PandasTools.LoadSDF(sdf)
    df = _calculate_descs(df,checks)
    return df

def read_smi(smi, checks):
    df = pd.read_csv(smi, delimiter="\t", names=['SMILES', 'ID'])
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol='SMILES')
    df = _calculate_descs(df,checks)
    return df
    