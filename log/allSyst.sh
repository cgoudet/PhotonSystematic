server=`pwd`
cd ${server} 
ulimit -S -s 100000 
LD_LIBRARY_PATH=/sps/atlas/c/cgoudet/Hgam/FrameWork/RootCoreBin/lib:/sps/atlas/c/cgoudet/Hgam/FrameWork/RootCoreBin/bin:$LD_LIBRARY_PATH 
cd /sps/atlas/c/cgoudet/Hgam/FrameWork/RootCoreBin/ 
source local_setup.sh 
cd ${server} 
cp -v /sps/atlas/c/cgoudet/Hgam/FrameWork/RootCoreBin/obj/x86_64-slc6-gcc49-opt/PhotonSystematic/bin/TestSyst . 
cp -v /sps/atlas/c/cgoudet/Hgam/Inputs/MxAOD_h013_Full/ntuple/* .
TestSyst --mode 1 --inConfFile /sps/atlas/c/cgoudet/Hgam/FrameWork/PhotonSystematic/data/FitFull.boost --outFileName allSyst/ group.phys-higgs.9486788._000006.MxAOD.root group.phys-higgs.9486788._000011.MxAOD.root group.phys-higgs.9486788._000016.MxAOD.root group.phys-higgs.9486788._000019.MxAOD.root group.phys-higgs.9486788._000023.MxAOD.root group.phys-higgs.9486788._000028.MxAOD.root group.phys-higgs.9486788._000029.MxAOD.root group.phys-higgs.9486788._000035.MxAOD.root group.phys-higgs.9486788._000045.MxAOD.root group.phys-higgs.9486815._000003.MxAOD.root group.phys-higgs.9486842._000007.MxAOD.root group.phys-higgs.9486842._000008.MxAOD.root group.phys-higgs.9486856._000004.MxAOD.root group.phys-higgs.9486881._000007.MxAOD.root group.phys-higgs.9486881._000011.MxAOD.root group.phys-higgs.9486881._000012.MxAOD.root group.phys-higgs.9486881._000018.MxAOD.root group.phys-higgs.9486881._000019.MxAOD.root group.phys-higgs.9486881._000020.MxAOD.root group.phys-higgs.9486881._000022.MxAOD.root group.phys-higgs.9486881._000028.MxAOD.root group.phys-higgs.9486881._000029.MxAOD.root group.phys-higgs.9486881._000030.MxAOD.root group.phys-higgs.9486881._000036.MxAOD.root group.phys-higgs.9486881._000041.MxAOD.root group.phys-higgs.9486881._000042.MxAOD.root group.phys-higgs.9486881._000049.MxAOD.root group.phys-higgs.9486881._000051.MxAOD.root group.phys-higgs.9486881._000056.MxAOD.root group.phys-higgs.9486881._000063.MxAOD.root group.phys-higgs.9486881._000066.MxAOD.root group.phys-higgs.9486881._000070.MxAOD.root
ls output/*
mkdir /sps/atlas/c/cgoudet/Hgam/FrameWork/Results/allSyst/
cp -v allSyst/* /sps/atlas/c/cgoudet/Hgam/FrameWork/Results/allSyst/.
