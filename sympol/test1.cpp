#include <iostream>
#include <fstream>

#include "polyhedron.h"
#include "polyhedronio.h"
#include "polyhedrondatastorage.h"
#include "automorphismcomputation.h"
#include "raycomputationlrs.h"
#include "recursionstrategy.h"
#include "revision.h"

using namespace sympol;
using namespace std;

int main(int argc, char* argv[]) {
  if (argc < 2)
    return 2;
  
  std::cout << "SymPol v0.0." << SVN_REVISION << " with lrs 4.2c and NTL 5.5.2" << std::endl;
  yal::ReportLevel::set(yal::DEBUG);
  
  std::string inputFile(argv[1]);
  
  ifstream myfile(inputFile.c_str());
  if (!myfile.is_open()) {
      return -1;
  }
  Polyhedron* poly = PolyhedronIO::read(myfile);
  myfile.close();
  
  RayComputationLRS* lrs = new RayComputationLRS();
  lrs->initialize();
  
  if (poly) {
    cout << "Poly" << endl;
    Polyhedron::RowIterator itBegin = poly->rowsBegin();
    Polyhedron::RowIterator itEnd = poly->rowsEnd();
    
    for (Polyhedron::RowIterator it = itBegin; it != itEnd; ++it) {
      cout << (*it).index() << " -- " << (*it) << endl;
    }
    
    boost::shared_ptr<PermutationGroup> autom(AutomorphismComputation::computeRestrictedIsomorphisms(*poly));
    if (autom)
      cout << "BSGS : " << *autom;
    else
      cerr << "no autom" << endl;
    
    /*
    std::vector<RayData> rays;
    lrs->dualDescription(*poly, rays);
    
    BOOST_FOREACH(const RayData& r, rays) {
      cout << " ray   " << r.ray << endl;
    }
    */
    
    std::list<RayData> rays2;
    //RecursionStrategyDirect* rsd = new RecursionStrategyDirect();
    RecursionStrategy* rsd = new RecursionStrategyADM();
    rsd->enumerateRaysUpToSymmetry(lrs, *poly, *autom, rays2);
    delete rsd;
    
    std::cout << "V-representation" << std::endl;
    PolyhedronIO::write(rays2, poly->homogenized(), cout);
    
    /*
    cout << "  ~~~" << endl;
    Face f(4);
    f[0] = 1;
    f[1] = 1;
    Polyhedron p2 = poly->supportCone(f);
    Polyhedron::RowIterator itBegin2 = p2.rowsBegin();
    Polyhedron::RowIterator itEnd2 = p2.rowsEnd();
    for (Polyhedron::RowIterator it = itBegin2; it != itEnd2; ++it) {
      cout << (*it).index() << " -- " << (*it) << endl;
    }*/
  } else {
    cerr << "could not read poly\n";
  }
 
  lrs->finish();
  delete lrs;
  
  
  PolyhedronDataStorage::cleanupStorage();
  
  return 0;
}