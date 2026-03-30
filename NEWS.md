
# Rpdb 2.4.4

## Refactoring
- Renaming "conect" -> "connect";
  - Renamed class to "connect";
  - Renamed function arguments to "connect";
  - Renamed functions to "connect";

## New Functionality
- Function chains: extract chain-codes;
- Function select: Select groups of atoms;
- Function residues: hidden/undocumented;
- Function translate: shift observer;
- Function addBBox: tight bounding box;
- S3 method: distances.connect;

## Enhancement:
- PDB:
  - Properly import chemical symbols of elements (for newer PDB format);
  - Improved handling of atom names (older versions of PDB format);
  - Re-index ElementID (always);
- Function addPBCBox:
  - new argument scale: includes also relative shift;
  - new argument col: set colour of box;
  - new argument alpha: set transparency;
  - x can be also a PDB molecule;
- Visualization:
  - Improved: Protein backbone;
  - Improved: Nucleic acid backbone;
  - [fix] Nucleic acids: Detect both "X\*" and "X'" codes;
  - Improved: Brute force connectivity;


## TODO

- Extract all useful code and ideas from the discussion thread:
  - Root: https://stat.ethz.ch/pipermail/r-help/2023-October/478359.html
  - Example: https://stat.ethz.ch/pipermail/r-help/2023-October/478372.html

=====

# Rpdb 2.4.3

- Improved visualization:
  - Connect protein backbone: basic code;

## Refactoring:
- Rename: conect -> connect;
  - started some work;


# Rpdb 2.4.1

## New Functions

- dist.point: distances to a given point;
- centres.ppRoll: centres computed using a rolling window over a peptide chain;


## Updates to Existing Functions

- Function read.pdb:
  - Bug fixed: process only FieldName CRYST1;
  - New arg: verbose = TRUE;
  - Warn if pdb is marked as obsolete;
  - Multi-line titles: import only the text;
  - Extract resolution from pdb;
  - Renaming arg CRYST1 to CRYSTAL;
- Function is.crystal: as replacement to is.cryst1;
- Function visualize.coords: starting some refactoring / optimization;
- Function toSymbols: started refactoring;
- [started] Renaming function cryst1 to crystal;

## Deprecated

- Class cryst1: replaced by class crystal;


# Rpdb 2.3.4

- new maintainer;

## Updates of Existing Functions

- code updated to use rgl >= 1.1.3;

## Documentation

- corrected various bugs with Roxygen2;
- improved compliance with CRAN-checks;
