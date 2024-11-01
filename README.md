# TemSOMap
Mapping lineage-resolved scRNA-seq data with spatial transcriptomics using TemSOMap
<img width="961" alt="Screenshot 2024-10-30 152639" src="https://github.com/user-attachments/assets/270463ac-9688-43cf-a4fd-83c5eee7e83c">

**Table of contents**

* [Installation](#Installation)
* [Usage](#Usage)
* [Results](#results)
* [Contact](#contact)

## Installation
The Python package for TemSOMap can be installed via PyPI as:
`pip install TemSOMap`

#### Please Cite

```
Citation available soon.
```

## Usage
Example usages of LinRace can be refered in vignettes.
%1. [Running TemSOMap on an example, simulated dataset](test/test_temso.ipynb)

## Results

The following figure shows an example application of TemSOMap, integrating a Stereo-seq and a lineage tracing dataset for E9.5 mouse embryo:
<img width="506" alt="Screenshot 2024-10-31 202019" src="https://github.com/user-attachments/assets/7602034f-0fd8-48b4-80d9-596a22c15a14">


**a**. 2-D visualization of the input Stereo-seq mouse embryo data (spot level) and the TemSOMap-inferred mouse embryo data (single-cell level). Colors represent cell types from the Stereo-seq annotation.
**b**. Spatiotemporal analysis on the cell fate specifications on the cell lineage. Each embryo plot represents the spatial distribution of cells for each subtree, with pie plots showing the cell type percentages of the leaves.
**c**.  Comparing spatial maps of gene expressions of the reference ST data (left), inferred cell-level data (middle), and inferred spot-level data (right). Foxc1 and Foxa2 are used in the training of TemSOMap; Ankrd1 and Meox1 are masked from the training of TemSOMap and their spatial distribution is predicted using TemSOMap.

## Contact
Other datasets used in the paper can be requested and
GitHub issues are welcomed.
It is also possible to send email to the main author
`Xinhai Pan (xpan78 at gatech.edu)`.
