import argparse
import binascii

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def run(input, separator, header=False, fragments=False):
    formulas = {}
    count = 0

    if fragments:
        print("Generating formula for each fragment")
    else:
        print("Generating formula for entire molecule")

    with open(input, "rt") as file:
        if header:
            file.readline()
        for line in file:
            if count % 1000000 == 0:
                print('... processed', count, len(formulas))
            tokens = line.split(separator)
            smiles = tokens[0]
            ok = True
            if fragments:
                parts = sorted(smiles.split('.'))
                mol_formula = None
                for part in parts:
                    mol = Chem.MolFromSmiles(part)
                    if not mol:
                        ok = False
                        continue

                    f = rdMolDescriptors.CalcMolFormula(mol)
                    if mol_formula:
                        mol_formula += '.' + f
                    else:
                        mol_formula = f
            else:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    mol_formula = rdMolDescriptors.CalcMolFormula(mol)
                else:
                    ok = False

            if ok:
                count += 1
                #crc = binascii.crc32(str.encode(dot_formula))
                #print(dot_formula, crc)
                accumulate(formulas, mol_formula)


    if formulas:
        f_counts = {}
        for f, c in formulas.items():
            accumulate(f_counts, c)
        for k in sorted(list(f_counts.keys())):
            print(k, f_counts[k])
        print("Formula size:", len(formulas))

    print("Processed", count, "molecules")

def accumulate(dict, item):
    value = dict.get(item, 0)
    dict[item] = value + 1


def main():
    parser = argparse.ArgumentParser(description="Analyse molecules")

    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-s", "--separator", default="\t")
    parser.add_argument("-l", "--header-line", action="store_true")
    parser.add_argument("-f", "--fragments", action="store_true")

    args = parser.parse_args()

    run(args.input, args.separator, header=-args.header_line, fragments=args.fragments)


if __name__ == "__main__":
    main()