from rdkit import Chem
from rdkit.Chem import rdBase
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit import RDConfig
from PIL import Image
import os

#file1 = "N#Cc1ccccc1Cn3c(N2CCCCC2)cc(=O)[nH]c3=O"

def moltosvg(file1,molSize=(400,400),kekulize=True):
    mol = Chem.MolFromSmiles(file1)
    AllChem.Compute2DCoords(mol)
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    img = Draw.MolToImage(mc)
    return img