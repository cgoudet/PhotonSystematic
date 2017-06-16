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
If the dataset is splitted in parts with distinct NP, create a list container1(2,3, ...).txt for each subset.


### Create ntuple from MxAOD
- Call python PhotonSystematic/python/FillNtuple.py --directory <path_to-version_directory> --mode 1
mode 1 launch independently each dataset on the batch, mode 0 launch successively all the datasets locally.

In each directory, a file FillRecord.txt saves the name of the file already launched. 
If one wants to relauch a file, it must be removed from this text file.

As of 20170614, FillNtuple recognises if the set of NP is splitted.
If the filename contains PhotonAllSysI, listContainersI.txt is used.

For each version, one has to check that the right branch is used for weight and caetgorisation.
This can be changed in PhotonSystematic/Root/FillNtuple.cxx/FillEtry(L.180).

Finished files will appear in the directory from which the routine is launched.

### Copy ntuple to directory
- Call python PhotonSystematic/python/FillNtuple.py --mode 2
- You can remove the FillNtuple directories in the current directory afterward

### Run the fitting procedure
- Call TestSyst <rootFile> --inConfFile <configurationFile> --mode 1 --outFileName <outDirectory>

The directory <outDirectory> is created if it doesn't exist. All the output files will be printed in this directory.

Configuration file must be in boost format. The available options are in FitSystematic::Configure
