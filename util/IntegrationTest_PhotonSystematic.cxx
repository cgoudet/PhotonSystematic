#include "PhotonSystematic/FitTree.h"
#include <stdexcept>
using std::runtime_error;
using std::cout;
using std::endl;

int main() {

  bool TestFitTree();
  if ( !TestFitTree() ) throw runtime_error( "Failed : TestFitTree" );
  cout << "All Tests OK" << endl;
  return 0;
}
