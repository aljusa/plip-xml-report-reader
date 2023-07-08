import xml.etree.ElementTree as ET
import pandas as pd
from collections import Counter
import sys

aa = {'ASP': ['acidic', 'acyclic', 'charged', 'medium', 'negative', 'polar', 'surface'],
        'GLU': ['acidic', 'acyclic', 'charged', 'large', 'negative', 'polar', 'surface'],
        'ALA': ['acyclic', 'aliphatic', 'buried', 'hydrophobic', 'neutral', 'small'],
        'GLY': ['acyclic', 'aliphatic', 'neutral', 'small', 'surface'],
        'LEU': ['acyclic', 'aliphatic', 'buried', 'hydrophobic', 'large', 'neutral'],
        'SER': ['acyclic', 'neutral', 'polar', 'small', 'surface'],
        'VAL': ['acyclic', 'aliphatic', 'buried', 'hydrophobic', 'medium', 'neutral'],
        'THR': ['acyclic', 'medium', 'neutral', 'polar', 'surface'],
        'LYS': ['acyclic', 'basic', 'charged', 'large', 'positive', 'surface'],
        'ILE': ['acyclic', 'aliphatic', 'buried', 'hydrophobic', 'large', 'neutral', 'polar', 'surface'],
        'ASN': ['acyclic', 'medium', 'neutral', 'polar', 'surface'],
        'GLN': ['acyclic', 'large', 'neutral', 'polar', 'surface'],
        'CYS': ['acyclic', 'buried', 'medium', 'neutral', 'polar'],
        'MET': ['acyclic', 'buried', 'hydrophobic', 'large', 'neutral'],
        'HIS': ['aromatic', 'basic', 'charged', 'cyclic', 'large', 'neutral', 'polar', 'positive', 'surface'],
        'PHE': ['aromatic', 'buried', 'cyclic', 'hydrophobic', 'large', 'neutral'],
        'TRP': ['aromatic', 'buried', 'cyclic', 'hydrophobic', 'large', 'neutral'],
        'TYR': ['aromatic', 'cyclic', 'large', 'neutral', 'surface'],
        'ARG': ['basic', 'charged', 'large', 'polar', 'positive', 'surface'],
        'PRO': ['cyclic', 'hydrophobic', 'medium', 'neutral', 'surface']}

ii = {'hydrophobic_interaction': 'hphobic',
        'hydrogen_bond':'hbond',
        'water_bridge':'water',
        'salt_bridge':'sbridge',
        'pi_stack':'pistack',
        'pi_cation_interaction':'pication',
        'halogen_bond':'halogen',
        'metal_complex':'metal'}

ll = {'num_heavy_atoms': 'heavy',
        'num_hbd':'hbd',
        'num_unpaired_hbd':'unpairedhbd',
        'num_unpaired_hba':'unpairedhba',
        'num_hba':'hba',
        'num_hal':'hal',
        'num_unpaired_hal':'unpairedhal',
        'num_aromatic_rings':'rings',
        'num_rotatable_bonds':'rotatable',
        'molweight':'molweight',
        'logp':'logp'}

def get_bsresidues(bs_residues):
    bsresidues = Counter()
    for element in bs_residues:
        if element.attrib['contact'] == "True":
            bsresidues.update(aa[element.attrib['aa']])
    return dict(bsresidues)

def get_ligproperties(lig_properties):
    ligproperties = {}
    for element in lig_properties:
        ligproperties[ll[element.tag]] = element.text

    return ligproperties

def get_interactions(_interactions):
    interactions = Counter()
    for element in _interactions:
        for subelement in element:
            intertype = ii[subelement.tag]
            resnr = subelement.findall('resnr')[0].text
            restype = subelement.findall('restype')[0].text
            reschain = subelement.findall('reschain')[0].text
            interactions.update([intertype + "_" + restype.title() + resnr + reschain])

    return dict(interactions)

def parse_bindingsite(bindingsite):
    for element in bindingsite:
        if element.tag == 'lig_properties':
            ligproperties = get_ligproperties(element)
        elif element.tag == 'bs_residues':
            bsresidues = get_bsresidues(element)
        elif element.tag == 'interactions':
            interactions = get_interactions(element)

    return ligproperties, bsresidues, interactions

def parse_report(xml):
    tree = ET.parse(xml)
    root = tree.getroot()
    ligands = []
    for i in root:
        if i.tag == 'bindingsite':      
          name = i.find("identifiers/members/member").text
          ligproperties, bsresidues, interactions  = parse_bindingsite(i)
          ligands.append((name, {'Receptor': xml[:-4],**ligproperties, **bsresidues, **interactions}))

    return ligands


def main():
    report_list = sys.argv[1] # reports list
    outfile = sys.argv[2] # output file

    ligands_all = []
    for file in open(report_list, "r").readlines():
        report_file = file.strip()
        print("Parsing " + report_file)
        ligands = parse_report(report_file)
        ligands.sort(key = lambda x: x[0].split(":")[1:])
        ligands_all.append(ligands)
   
    report = {}    
    for ligands in ligands_all:
        for name, ligand in ligands:
            report[name] = ligand
    df = pd.DataFrame.from_dict(report)
    df = df.transpose()
    df = df.rename_axis('Ligand')
    df.to_csv(outfile)

main()
