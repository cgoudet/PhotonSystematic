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

Each file is usually about 4GB heavy. 
The full PhotonAllSys h015d production weights 415GB.

### Extract the name of all the EGamma NP in the model.
- Open a MxAOD in root (do not do any rcSetup in the session) and CollectionTree->Print(); > ../<version>_<EGamModel>.txt
- Call PhotonSystematic/python/FindBranches.py <path_to version directory>. 
This will generate alistContainers.txt files with all the branches in the current file.
It is assumed that the branches of interest in the files are all the same.
If the dataset is splitted in parts with distinct NP, the procedure should be repeated with of each set and the listconstainers files must be renamed listcontainer1(2,3, ...).txt.


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

A typical files contains 100k events.
This steps is usually run at 15Hz interactively (--mode 0) for a running time of about two hours.
On the batch, the code usually runs at lower rate (5-10Hz) and lasts for 3 to four hours.
Given the large number of files, it is usually run of the batch (--mode 1).

### Copy ntuple to directory
- Add the file version to the python code responsible for copying. In PhotonSystemati/python/FillNtuple.py, in the function MoveSample, you can add a keyword which define the version.
- Call python PhotonSystematic/python/FillNtuple.py --mode 2
- You can remove the FillNtuple directories in the current directory afterward

### Run the fitting procedure
- Call TestSyst 'rootFile' --inConfFile 'configurationFile' --mode 1 --outFileName 'outDirectory'

The directory <outDirectory> is created if it doesn't exist. All the output files will be printed in this directory.

Configuration file must be in boost format. The available options are in FitSystematic::Configure

For h015d production, the total fitting time is estimated at 35min : 30mins to read ntuples and fill datasets, 5min for the fits.

### Create the plots
- python PhotonSystematic/python/PlotTests.py 'resultFile'

By default, plots for the full set of systematics are created.
If the name of the directory contains 'Merge', the postmerge results will be used.
Carefull, the TestSyst program use the directory name as prefix for all the files and this convention is used by PlotTest.
If one wants just to rename a directory, mass renaming all the files inside is also necessary.
