#ifndef PhotonSystematic_TutorialAnalysis_H
#define PhotonSystematic_TutorialAnalysis_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include <map>
#include <string>
using std::string;
using std::map;

class TutorialAnalysis : public HgammaAnalysis
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
private:
  // Tree *myTree; //!
  // TH1 *myHist; //!


  TTree *m_output_tree; //!
  map<string, double> m_mapTree;

public:
  // this is a standard constructor
  TutorialAnalysis() { }
  TutorialAnalysis(const char *name);
  virtual ~TutorialAnalysis();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();


  // this is needed to distribute the algorithm to the workers
  ClassDef(TutorialAnalysis, 1);
};

#endif // PhotonSystematic_TutorialAnalysis_H
