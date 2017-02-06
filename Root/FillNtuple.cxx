#include "PhotonSystematic/FillNtuple.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/HGamVariables.h"
#include "PlotFunctions/Foncteurs.h"
#include "PlotFunctions/SideFunctions.h"
using namespace ChrisLib;
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <algorithm>
using std::copy;
using std::for_each;
using std::set_difference;
using std::transform;
using std::all_of;
#include <exception>
using std::ostream_iterator;
using std::back_inserter;
using std::runtime_error;
#include <iterator>
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>

using std::set;
using std::string;
using std::vector;
using std::map;
using std::list;
// this is needed to distribute the algorithm to the workers
ClassImp(FillNtuple)



FillNtuple::FillNtuple(const char *name)
: HgammaAnalysis(name)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



FillNtuple::~FillNtuple()
{
  // Here you delete any memory you allocated during your analysis.
}



EL::StatusCode FillNtuple::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.

  //  histoStore()->createTH1F("m_yy", 60, 110, 140);

  m_containersName.insert( "" );
  m_debug=1;

  TFile *file = wk()->getOutputFile("MxAOD");
  m_outTree = new TTree("output","output");
  m_outTree->SetDirectory(file);

  string containerConfigName = config()->getStr("containerConfig").Data();
  DefineContainers( containerConfigName );
  LinkOutputTree();
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode FillNtuple::execute()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  HgammaAnalysis::execute();


  for ( auto itVar : m_analysisVariables ) m_branchLinks[itVar] = -99;

  int mcID=eventInfo()->mcChannelNumber(); //341000=ggH
  m_branchLinks["weight"] = lumiXsecWeight(10.,mcID,1) * var::weightCatCoup_Moriond2017();
  
  for (auto sys:getSystematics()) {
    if ( m_containersName.find(sys.name()) == m_containersName.end() ) continue;
    CP_CHECK("execute()",applySystematicVariation(sys));
    FillEntry( sys.name() );
  } 
  m_outTree->Fill();
  
  if ( m_debug==1 ) m_debug=0;
  return EL::StatusCode::SUCCESS;
}

//==============================================================
double FillNtuple::ReweightPtGgh( double const initPt ) {
  	if (initPt<20) return 1.11;
	if (initPt<45) return 1.11 - (initPt-20)/25*0.2; // -> 0.91
  	if (initPt<135) return 0.91 - (initPt-45)/90*0.36; // -> 0.55
  	return 0.55;
}

//==============================================================
void FillNtuple::DefineContainers( const std::string &containerConfig ) {

  vector<string> containersName;
    po::options_description desc("LikelihoodProfiel Usage");
    desc.add_options()
      ( "help", "Display this help message")
      ( "containerName", po::value<vector<string>>( &containersName )->multitoken(), "Names of the container to copy" )
      ( "commonVarName", po::value<vector<string>>( &m_commonVarsName )->multitoken(), "Names of the variables in containers to copy" )
      ;
    // create a map vm that contains options and all arguments of options       
    po::variables_map vm;
    std::ifstream ifs( containerConfig, std::ifstream::in );
    po::store(po::parse_config_file(ifs, desc), vm);
    po::notify(vm);

    for ( auto itName : containersName ) m_containersName.insert(itName);

    sort( m_commonVarsName.begin(), m_commonVarsName.end() );
    for ( auto vVar : m_commonVarsName ) {
      auto it = std::find( m_analysisVariables.begin(), m_analysisVariables.end(), vVar );
      if ( it==m_analysisVariables.end() ) throw runtime_error( "FillNtuple::DefineContainers : " + vVar + "is not a known variable." );
    }
    
    set_difference( m_analysisVariables.begin(), m_analysisVariables.end(), m_commonVarsName.begin(), m_commonVarsName.end(), back_inserter(m_systVarsName) );
}

//==============================================================
void FillNtuple::LinkOutputTree() {

  list<list<string>> inCombineNames( 2, list<string>() );
  copy( m_containersName.begin(), m_containersName.end(), back_inserter( *inCombineNames.begin() ) );
  copy( m_systVarsName.begin(), m_systVarsName.end(), back_inserter( inCombineNames.back() ) );

  list<string> outBranchesName;
  CombineNames( inCombineNames, outBranchesName );
  outBranchesName.insert( outBranchesName.end(), m_commonVarsName.begin(), m_commonVarsName.end() );

  for ( auto it = outBranchesName.begin(); it!=outBranchesName.end(); ++it ) 
    m_outTree->Branch( it->c_str(), &m_branchLinks[*it] );


}

//==============================================================
void FillNtuple::FillLink( const string &inName, const string &outName, const map<string,double> &vars ) {
  map<string,double>::iterator citVarVal = m_branchLinks.find(outName);
  map<string,double>::const_iterator citLocalVal = vars.find(inName);
  if ( m_debug==1 ) {
    assert( citVarVal!=m_branchLinks.end() );
    assert( citLocalVal!=vars.end() );
  }
  if ( inName!=outName 
       || citVarVal->second == -99 
       ) citVarVal->second = citLocalVal->second;

}

//==============================================================    
bool FillNtuple::FillEntry( const string &systName ) {
  map<string,double> vars;
  double dummy = var::m_yy();
  vars["m_yy"] = dummy==-99 ? dummy : dummy/1e3;
  vars["weight"] = var::weight();
  double dummyVar = var::catCoup_Moriond2017();
  vars["catCoup"] = dummyVar==-1 ? -99 : dummyVar;

  for ( const string itSyst:m_systVarsName ) {
    const string outName { ( systName =="" ? "" : systName+"_") +itSyst };
    FillLink( itSyst, outName, vars );
  }
  for ( const string itCommon:m_commonVarsName ) FillLink( itCommon, itCommon, vars );
  return true;
}
