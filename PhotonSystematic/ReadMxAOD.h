#ifndef READMXAOD_H
#define READMXAOD_H

#include "TFile.h"
#include "xAODEventInfo/EventInfo.h"
#include <vector>
#include <string>
#include <map>
#include <list>
using std::map;
using std::string;
using std::list;
using std::vector;
/*\brief Read an MxAOD of systematic variation a print a set a variable for each variation into a ntuple
  \param inConfFileName Input configuration file in boost format. (more below)
  \param debuf Debug mode. 3 modes are defined : 0 user mode, 1 medium printing over full MxAOD inputs, 2 maximum print over 10 events

  
 */
int ReadMxAOD( const string &inConfFileName, int debug = 0 );

/*\brief Reweight diphton events according to the pt of the system
  \param initPt Truth Pt of the diphoton system.
  Macro given by Dag. Tested.
 */
double ReweightPtGgh( double const initPt );

/*\brief Find in the MxAOD the histograms which contains the sum of weight to use on the events
  \param inFile in which to search for the histogram
  \return Name of the noDalitz_weighted histogram
  Tested
 */
string FindNoDalitzHist( const TFile *inFile );

/*\brief Fill a map with the read values of a tevent
 */
bool FillMapFromEventInfo( const string &outName, map<string,double> &mapVal, const xAOD::EventInfo* eventInfo, double commonWeight=1, bool isCommon=0 );

/*\brief Look in a map for variables that may have different values in different branches and remove them
  \param dumplicateListName list all the variables
  \param mapVal values
  \param Possible default values of the variables.

  The code assume that the map keys are in the form branch_variable with no underscore in variable except for _yy. 
  Tested
*/
void UpdateDuplicateList( list<string> &duplicateListName, const map<string, double> &mapVal, const map<string, double> &defaultValues );

/*\brief Sum the values of the 3rd bin in noDalitz_weighted histograms in all input files, defenciating process thourh name
  \param rootFilesName names of the input files
  \param mapVal values as a function of process

  The process of a file is obtained by looking for the process name in the file name.
  Tested.
 */
void TotalSumWeights( const vector<string> &rootFilesName, map<string,double> &datasetWeights );

/*\brief Extract variable name in a string
\param inName  string to parse
The code assumes the string writtent as branchName_varaible wit no underscore in variable execpt for final _yy.
Tested.
*/
string ExtractVariable( const string &inName );
#endif
