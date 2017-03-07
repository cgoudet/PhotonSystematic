#ifndef PhotonSystematic_FillNtuple_H
#define PhotonSystematic_FillNtuple_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"


#include <vector>
#include <string>
#include <map>
#include <set>

class FillNtuple : public HgammaAnalysis
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
private:
  //outputs
  TTree *m_outTree; //!
  int m_debug; //!

  std::set<std::string> m_containersName;//!
  std::vector<std::string> m_commonVarsName;//!
  std::vector<std::string> m_systVarsName;//!
  std::map<std::string,double> m_branchLinks;//!

  double ReweightPtGgh( double const initPt ); //!
  void DefineContainers( const std::string &containerConfig );
  void LinkOutputTree();//!
  bool FillEntry( const std::string &systName );//!
  void FillLink( const std::string &inName, const std::string &outName, const std::map<std::string,double> &vars ); //!

  const std::list<std::string> m_analysisVariables { "catCoup", "m_yy", "weight", "catCoupBDT"  }; //!
  const std::list<std::string> m_processes { "ggH", "VBF", "WH", "ZH", "ttH", "bbH125_yb2", "bbH125_ybyt", "tWH", "tHjb" };
  //  static const std::list<std::string> vect = { "DPhi_yy","catCoup","catXS","catXSPhi","m_yy", "pt_yy","weight","weightXS"};

public:
  // this is a standard constructor
  FillNtuple() { }
  FillNtuple(const char *name);
  virtual ~FillNtuple();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();


  // this is needed to distribute the algorithm to the workers
  ClassDef(FillNtuple, 1);
};




#endif // PhotonSystematic_FillNtuple_H
