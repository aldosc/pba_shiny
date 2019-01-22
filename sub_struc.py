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
from rdkit.Chem import rdFMCS
from collections import namedtuple
from cairosvg import svg2png

def make_subst(file1,file2):
    res_df = pd.DataFrame(file1)
    list_smiles1=list(file1.canonical_smiles)
    list_smiles1.append(file2)
    list_smiles=[Chem.MolFromSmiles(x) for x in list_smiles1]
    list_ids=list(res_df.name)
    list_ids.append("My_query")
    
    mhs = list_smiles
    [mhs[j].SetProp('_Name', str(list_ids[j])) for j in range(0, len(mhs))];
    mcs = rdFMCS.FindMCS(mhs,threshold=0.8,completeRingsOnly=True,ringMatchesRingOnly=True)
    patt = Chem.MolFromSmarts(mcs.smartsString)

    AllChem.Compute2DCoords(patt)


    matchingMols = [x for x in mhs if x.HasSubstructMatch(patt)]
    legends=[x.GetProp("_Name") for x in matchingMols]

    for m in matchingMols: AllChem.GenerateDepictionMatching2DStructure(m,patt)

    ##code for visualizing the sub-structure
    hats = []
    hbnds = []
    for mm in matchingMols:
        ats = mm.GetSubstructMatch(patt)
        hats.append(ats)
        bnds = []
        for bnd in mm.GetBonds():
            if bnd.GetBeginAtomIdx() in ats and bnd.GetEndAtomIdx() in ats:
                bnds.append(bnd.GetIdx())
        hbnds.append(bnds)

    molsPerRow = 3
    nRows = len(matchingMols)//molsPerRow
    if len(matchingMols)%molsPerRow:
        nRows+=1
    panelx = 250
    panely = 260
    canvasx = panelx * molsPerRow
    canvasy = panely * nRows
    drawer = rdMolDraw2D.MolDraw2DSVG(canvasx,canvasy,panelx,panely)
    drawer.DrawMolecules(matchingMols,highlightAtoms=hats,highlightBonds=hbnds, legends=legends)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    SVG(svg.replace("svg:",""))
    output = svg2png(bytestring=svg)
    return svg
