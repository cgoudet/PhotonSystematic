\mainpage

# PhotonSystematic

## Steps to produce signal shape systematics using the HGam MxAOD. 
### Download the Photon(All)Sys datasets into a folder <prod>_<EGamModel>/MxAOD.
The datasets which must be downloaded are : 
- mc15c.MGPy8_tHjb125_yt_plus1
- mc15c.PowhegPy8_NNLOPS_ggH125
- mc15c.PowhegPy8_NNPDF30_VBFH125
- mc15c.PowhegPy8_WmH125J
- mc15c.PowhegPy8_WpH125J
- mc15c.PowhegPy8_ZH125J
- mc15c.PowhegPy8_ggZH125
- mc15c.aMCnloHwpp_tWH125_yt_plus1.
- mc15c.aMCnloPy8_bbH125_yb2
- mc15c.aMCnloPy8_bbH125_ybyt
- mc15c.aMCnloPy8_ttH125


### Extract the name of all the EGamma NP in the model.
- Open a MxAOD in root (do not do any rcSetup in the session) and CollectionTree->Print(); > ../<version>_<EGamModel>.txt
- Call PhotonSystematic/python/FindBranches.py <path_to version directory>. 
This will generate the listContainers.txt files with all the branches in the current file.
The branches should be harmonized between all input files