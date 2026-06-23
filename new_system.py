#!/usr/bin/env python3
"""Scaffold a YAML systems block for a new hub gene and its partners.

  python new_system.py --hub MYHUB --partner PARTNERA --partner PARTNERB

Prints a ready-to-fill YAML block and lists the structure files you must provide.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
from mavis_v7.config_io import scaffold_system_yaml


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hub", required=True, help="Hub gene (the one your variants are in)")
    ap.add_argument("--partner", action="append", required=True, dest="partners",
                    help="A binding partner (repeat for several)")
    args = ap.parse_args()
    text, files = scaffold_system_yaml(args.hub, args.partners)
    print(text)
    print("# --- structure files to place under your --structures directory ---")
    for f in files:
        print(f"#   {f}")
    print("#   (+ a .cif alongside each AlphaFold .pdb if your PDB B-factors are zeroed)")
    print("#\n# Then run:")
    print("#   python run.py --config your_systems.yaml --variants your_variants.csv "
          "--structures <dir> --out results/run")


if __name__ == "__main__":
    main()
