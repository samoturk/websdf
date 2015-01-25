from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools

def read_sdf(sdf):
    df = PandasTools.LoadSDF(sdf)
    return df
    