#ifndef FITTREE_H
#define FITTREE_H

#include "PlotFunctions/MapBranches.h"

#include "RooDataSet.h"
#include "RooRealVar.h"

#include <string>
#include <map>
#include <list>

int FitTree( std::string inConfFileName );
void FillInitialValuesFitParam( std::map<std::string,std::vector<double>> &mapInitValues );
void FillEntryDataset( const std::list<std::string> &NPName, 
		       const MapBranches &mapBranch, 
		       map<std::string,vector<RooDataSet*>> &mapSet,
		       map<std::string,RooRealVar> &observables,
		       const std::string &catVar );
#endif
