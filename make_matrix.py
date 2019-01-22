import numpy as np
import pandas as pd

from lxml import etree

from IPython.display import Image, SVG

from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit import rdBase

from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

def make_matrix(file1, file2):
    res_df = pd.DataFrame(file1)
    list_smiles1=list(file1.canonical_smiles)
    list_smiles1.append(file2)
    list_smiles=[Chem.MolFromSmiles(x) for x in list_smiles1]
    list_ids=list(res_df.name)
    list_ids.append("My_query")
    my_fps = [FingerprintMols.FingerprintMol(x) for x in list_smiles]

    dists = []
    simil = []
    nfps = len(my_fps)
    for j in range(0,nfps):
        simil.append(DataStructs.BulkTanimotoSimilarity(my_fps[j],my_fps))
        res_dis = DataStructs.BulkTanimotoSimilarity(my_fps[j],my_fps,returnDistance=1)
        dists.append([1-x for x in res_dis])
    
    simil_mat = np.array(simil)
    dist_mat = np.array(dists)
    df_dist = pd.DataFrame(dist_mat)
    df_simil = pd.DataFrame(simil_mat)
    df_simil.columns = list_ids
    df_simil.index = list_ids
    return df_simil
