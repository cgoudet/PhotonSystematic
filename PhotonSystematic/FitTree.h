#ifndef FITTREE_H
#define FITTREE_H

#include <string>
#include <map>

int FitTree( std::string inConfFileName );
void FillInitialValuesFitParam( std::map<std::string,std::vector<double>> &mapInitValues );

#endif
