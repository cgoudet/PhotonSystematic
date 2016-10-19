#ifndef READMXAOD_H
#define READMXAOD_H

#include <string>
using std::string;
#include "TFile.h"

int ReadMxAOD( string inConfFileName, int debug = 0 );
double ReweightPtGgh( double initPt );
string FindNoDalitzHist( TFile *inFile );
#endif
