from __future__ import print_function
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from cairosvg import svg2png

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
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    SVG(svg.replace("svg:",""))
    svg2png(bytestring=svg,write_to='molecule2_out_pba.png')
    return svg