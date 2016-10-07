#include "PhotonSystematic/TutorialAnalysis.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/HGamVariables.h"

#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector;

// this is needed to distribute the algorithm to the workers
ClassImp(TutorialAnalysis)



TutorialAnalysis::TutorialAnalysis(const char *name)
: HgammaAnalysis(name)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



TutorialAnalysis::~TutorialAnalysis()
{
  // Here you delete any memory you allocated during your analysis.
}



EL::StatusCode TutorialAnalysis::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.

  //histoStore()->createTH1F("m_yy", 60, 110, 140);

  TFile *file = wk()->getOutputFile("MxAOD");
  m_output_tree = new TTree("output","output");
  m_output_tree->SetDirectory(file);



  vector<string> var = { "mass", "coupCateg" };
  for ( auto vMap : var ) {
    m_mapTree[vMap] = 0;
    m_output_tree->Branch( vMap.c_str(), &m_mapTree[vMap] );
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TutorialAnalysis::execute()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  HgammaAnalysis::execute();

  xAOD::PhotonContainer photons = photonHandler()->getCorrectedContainer();
  xAOD::PhotonContainer photonsSel  = photonHandler()->applySelection(photons);
  if (photonsSel.size() < 2) return EL::StatusCode::SUCCESS;
  //  TLorentzVector h = photonsSel[0]->p4() + photonsSel[1]->p4();

  for ( auto iPhoton = 0; iPhoton<2; iPhoton++ ) {
    string sPhoton = string( TString::Format( "_%d", iPhoton ) );
    m_mapTree["eta" + sPhoton] = photonsSel[iPhoton]->eta();
    m_mapTree["phi" + sPhoton] = photonsSel[iPhoton]->phi();
    m_mapTree["pt" + sPhoton] = photonsSel[iPhoton]->pt();
  }

  // for (auto sys: getSystematics()) {
  //   if (sys.name() == "") continue;
  //   TString suffix = sys.name();
  //   applySystematicVariation(sys);
  //   xAOD::PhotonContainer sysPhotons  = photonHandler()->getCorrectedContainer();
  //   xAOD::PhotonContainer allSysPhotons  = photonHandler()->applySelection(sysPhotons);
    
  //   if (allSysPhotons.size() < 2) continue;
  //   TLorentzVector yy = allSysPhotons[0]->p4() + allSysPhotons[1]->p4();
  //   yy *= HG::invGeV;
  //   for ( auto iPhoton = 0; iPhoton<2; iPhoton++ ) {
  //     string sPhoton = string( TString::Format( "_%d", iPhoton ) );
  //     m_mapTree["pt" + sPhoton + "_" + string(suffix)] = allSysPhotons[0]->pt();
  //   }
  //   cout << "Systematic name : " << suffix << " "  << allSysPhotons[0]->pt() << endl;
  // }



  // for ( auto vMap : m_mapTree ) {
  //   cout << vMap.first << " " << vMap.second << endl;
  // }
  //  if ( m_output_tree->GetEntries() == 3 ) exit(0);
  m_output_tree->Fill();

  return EL::StatusCode::SUCCESS;
}
