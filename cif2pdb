import argparse
import warnings
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="folder with .cif files")
    parser.add_argument("-o", "--output", required=True, help="folder with .pdb files")
    args = parser.parse_args()

    in_dir = Path(args.input)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    cif_files = list(in_dir.rglob("*.cif"))
    parser_cif = MMCIFParser(QUIET=True)
    io = PDBIO()

    print(f"{len(cif_files)} .cif found. Converting...")

    for cif in cif_files:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", PDBConstructionWarning)
                structure = parser_cif.get_structure(cif.stem, cif)
            
            pdb_path = out_dir / f"{cif.stem}.pdb"
            io.set_structure(structure)
            io.save(str(pdb_path))
        except Exception as e:
            print(f"Error processing {cif.name}: {e}")

    print(f"Conversion ended. PDB files in: {out_dir}")

if __name__ == "__main__":
    main()
