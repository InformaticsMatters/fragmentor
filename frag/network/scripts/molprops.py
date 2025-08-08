import argparse
import gzip

from rdkit import Chem
from rdkit.Chem import ChemicalFeatures, rdMolDescriptors, Mol, rdmolops
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory

def run(input, output):
    with gzip.open(input, 'rt') as f:
        with gzip.open(output, 'wt') as out:
            count = 0
            for line in f:
                count += 1
                if (count % 100000) == 0:
                    print('Processed ', count)
                tokens1 = line.strip().split(',')
                synthon, core, fp = calculate(tokens1)
                tokens1.append(synthon)
                tokens1.append(core)
                tokens1.append(fp)
                out.write(','.join(tokens1) + '\n')


def calculate(tokens1):
    smiles1 = tokens1[0]
    label = tokens1[2]
    tokens2 = label.split('|')

    synthon = tokens2[1]
    core = tokens2[4]
    if synthon.count('Xe') == 1:
        mol = Chem.MolFromSmiles(smiles1)
        fp = calc_pharmfp(mol)
    else:
        fp = ""
    return synthon, core, fp


featFactory = ChemicalFeatures.BuildFeatureFactory('FeatureswAliphaticXenon.fdef')
sigFactory = SigFactory(featFactory, maxPointCount=2)
sigFactory.SetBins([(0, 2), (2, 5), (5, 8)])
sigFactory.Init()
sigFactory.GetSigSize()

def calc_pharmfp(mol):
    fp = Generate.Gen2DFingerprint(mol, sigFactory)
    fp = list(fp)
    return ";".join(map(str, fp))



def main():
    # Read in a edges.csv.gz and calculate additional properties
    parser = argparse.ArgumentParser(
        description="Generate edge properties"
    )
    parser.add_argument("-i", "--input", required=True, help="file to process")
    parser.add_argument("-o", "--output", required=True, help="output file")


    args = parser.parse_args()

    run(args.input, args.output)



if __name__ == "__main__":
    main()
