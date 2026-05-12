# Rpdb: Read, Write, Visualize and Manipulate PDB Files

- Reanimated the Rpdb package from the CRAN archive;
- Original author: Julien Idé;

## Changes

- see [News.md] for all changes;

## Usage

```
library(Rpdb)
# or: devtools::load_all(_path_)

# Download some legacy PDB files:
x = read.pdb(file.choose())
visualize(x, type="l")
```

**Example:**
```
# Download legacy PDB files from:
# https://www.rcsb.org/
x = read.pdb("7K3G.pdb")

# NMR Structure:
# Name of function drop.h: not yet stable;
x0H = drop.h(x)
con = connect(x0H) # Brute-force connectivity;
# Visualisation:
visualize(atoms(x0H), connect = con, type="l", pbc = F); addBBox(x)

# Labels:
res = residues(x0H, "A")
cols = rep("#84A232", nrow(res))
cols[res$resname %in% c("SER", "THR", "ASN")] = "red"
mark.pdb(x0H, ch = "A", col = cols, lwd = 2)

```

Note:
- Visualisation is incomplete;
- Started redesigning the Bounding box vs PBC box;
- Some functions are not yet visible in Rpdb, but can be accessed using: Rpdb:::__function__(...);

<!-- badges: start -->
[![R-CMD-check](https://github.com/discoleo/Rpdb/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/discoleo/Rpdb/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
