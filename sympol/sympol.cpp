// ---------------------------------------------------------------------------
//
// This file is part of SymPol
//
// Copyright (C) 2006-2010  Thomas Rehn <thomas@carmen76.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
// ---------------------------------------------------------------------------

#include "configuration.h"
#include "polyhedron.h"
#include "polyhedronio.h"
#include "raycomputationlrs.h"
#include "raycomputationcdd.h"
#include "automorphismcomputation.h"
#include "recursionstrategy.h"
#include "recursionstrategyidmadm.h"
#include "recursionstrategyidmadmlevel.h"
#include "recursionstrategyadmidmlevel.h"
#include "symmetrycomputationdirect.h"
#include "symmetrycomputationadm.h"
#include "symmetrycomputationidm.h"
#include "yal/logger.h"
#include "yal/usagestats.h"
#include "config.h"

// PermLib
#include <permlib/version.h>
#include <permlib/construct/schreier_sims_construction.h>

#include <csignal>
#include <string>
#include <sstream>
#include <fstream>
#include <list>
#include <boost/shared_ptr.hpp>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/erase.hpp>

namespace po = boost::program_options;
using namespace permlib;
using namespace std;
using namespace sympol;
using namespace yal;

#ifdef SIGUSR1
void stats(int param) {
  std::cout << "#ADM " << SymmetryComputationADM::instanceNumber() << std::endl;
}
#endif //SIGUSR1

void print_group(const boost::shared_ptr<PermutationGroup>& pg) {
	cout << pg->S.size() << endl;
	BOOST_FOREACH( const PERMptr &p, pg->S ) {
		ostringstream ssPerm;
		ssPerm << *p;
		string sPerm(ssPerm.str());
		boost::algorithm::replace_all(sPerm, ",", " ");
		boost::algorithm::replace_all(sPerm, ")(", ",");
		boost::algorithm::erase_all(sPerm, "(");
		boost::algorithm::erase_all(sPerm, ")");
		cout << " " << sPerm << endl;
	}
	cout << pg->B.size() << endl;
	BOOST_FOREACH( ulong beta, pg->B ) {
		cout << " " << (beta+1);
	}
	cout << endl;

}

uint correct_id(uint id, uint apexIndex) {
	return (id > apexIndex) ? (id - 1) : id;
}

int main(int argc, char* argv[]) {
	std::cout << "SymPol v0.1.1 and PermLib " << PERMLIB_VERSION << " with lrs 4.2c and cddlib 0.94f";
#if HAVE_NAUTY && HAVE_NTL
	std::cout << " and nauty and NTL";
#endif
	std::cout << std::endl;
	
	std::cout << "called as";
	for (int i = 0; i < argc; ++i) {
		std::cout << " " << argv[i];
	}
	std::cout << std::endl;
	
	LoggerPtr logger = Logger::getLogger("SymPol    ");
	
#ifdef SIGUSR1
	if (signal (SIGUSR1, stats) == SIG_IGN)
		signal (SIGUSR1, SIG_IGN);
#endif //SIGUSR1

	// declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("verbose,v", po::value<uint>()->default_value(0)->implicit_value(1), "increase verbosity")
		("time,t", "enable time measurement")
		("input-file,i", po::value<string>(), "input file  (.ine or .ext)")
		("automorphisms-only", "compute only restricted automorphisms and then exit")
		("no-automorphisms", "ignores given symmetry group and does not attempt to compute restricted automorphisms")
#if HAVE_NAUTY && HAVE_NTL
		("nauty", "use nauty graph isomorphism tool and NTL to compute restricted isomorphisms")
#endif
		("estimation-only,e", "compute only LRS estimation and then exit")
		("direct,d", "compute dual description directly")
		// removed defunct lattice code in rev.1207; needs to be rewritten
		//("lattice,l", "compute face lattice")
		("adm,a", po::value<double>(), "use adjacency decomposition method up to given estimate threshold")
		("adm-dim", po::value<uint>(), "use adjacency decomposition method up to given dimension threshold")
		("adm-incidence", po::value<uint>(), "use adjacency decomposition method up to given incidence number threshold")
		("idm-adm", po::value<vector<double> >()->multitoken(), "combined IDM,ADM strategy, expects two parameters: thresholdIDM thresholdADM")
		("idm-adm-level", po::value<vector<uint> >()->multitoken(), "combined IDM,ADM strategy, expects two parameters: levelIDM levelADM")
		("adm-idm-level", po::value<vector<uint> >()->multitoken(), "combined ADM,IDM strategy, expects two parameters: levelADM levelIDM")
		("cdd", "use cdd for core dual description conversion (EXPERIMENTAL)")
		("adjacencies", "records facet adjacencies (requires ADM method at level 0)")

		// serialization removed in rev.1185-1186 because it didn't work correctly
		//("dump", po::value<string>()->implicit_value(defaultDumpFilename), "create regularly dump file with which an interupted computation can be easily resumed")
		//("resume", po::value<string>()->implicit_value(defaultDumpFilename), "resume a previous interrupted computation")
	;
	
	po::options_description confOptions("Parameters");
	confOptions.add_options()
		("conf-lrs-estimates", po::value<uint>(&Configuration::getInstance().lrsEstimates)->default_value(Configuration::getInstance().lrsEstimates), "number of estimates that LRS performs")
		("conf-lrs-estimate-maxdepth", po::value<uint>(&Configuration::getInstance().lrsEstimateMaxDepth)->default_value(Configuration::getInstance().lrsEstimateMaxDepth), "maximal allowed depth of LRS estimates")
		("conf-compute-invariants", po::value<uint>(&Configuration::getInstance().computeInvariants)->default_value(Configuration::getInstance().computeInvariants), "maximal degree of permutation group invariant polynomials to use")
		("conf-compute-orbit-limit", po::value<uint>(&Configuration::getInstance().computeOrbitLimit)->default_value(Configuration::getInstance().computeOrbitLimit), "memory limit for computing full orbit (in megabytes)")
		("conf-intermediate-poly-fileprefix", po::value<string>(&Configuration::getInstance().intermediatePolyFilePrefix)->default_value(""), "prefix for filenames in which intermediate polyhedra are saved (useful for debugging and analyzing instances)")
	;
	desc.add(confOptions);

	po::positional_options_description pd;
	pd.add("input-file", -1);
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
			options(desc).positional(pd).run(), vm);
	po::notify(vm);
	
	if (vm.count("verbose"))
		yal::ReportLevel::set(static_cast<yal::LogLevel>(yal::WARNING + vm["verbose"].as<uint>()));
	
	if (vm.count("help")) {
		cout << desc << endl;
		return 1;
	}
	
	if (!vm.count("input-file")) {
		cout << desc << endl;
		return 1;
	}
	
	string inputFile = vm["input-file"].as<string>();
	
	ifstream myfile(inputFile.c_str());
	if (!myfile.is_open()) {
		YALLOG_ERROR(logger, inputFile << " is no valid input file");
		return -1;
	}
	std::list<boost::shared_ptr<PERM> > groupGenerators;
	std::vector<ulong> groupBase;
	Polyhedron* poly = PolyhedronIO::read(myfile, groupGenerators, groupBase);
	myfile.close();
	
	if (!poly) {
		YALLOG_ERROR(logger, "could not read polyhedron; exiting");
		return -1;
	}
	
	const bool bTimeMeasurement = vm.count("time") > 0;
	double timeVal = 0.0;
	if (bTimeMeasurement) {
		timeVal = UsageStats::processTimeUser();
	}
	
	// prepare computation instruments
	RayComputation *rayComp = 0;
	if (vm.count("cdd"))
		rayComp = new RayComputationCDD();
	else
		rayComp = new RayComputationLRS();
	
	rayComp->initialize();
	
	boost::shared_ptr<PermutationGroup> pg;
	if (!vm.count("no-automorphisms")) {
		if (groupGenerators.empty()) {
#if HAVE_NAUTY
			if (vm.count("nauty"))
				pg = AutomorphismComputation::computeRestrictedIsomorphismsNauty(*poly);
#endif
			if (!pg)
				pg = AutomorphismComputation::computeRestrictedIsomorphisms(*poly);
		} else {
			YALLOG_INFO(logger, "input file contains group generators; skip automorphism search");
			SchreierSimsConstruction<PERM, TRANSVERSAL> ssc(poly->realRowNumber());
			pg = boost::shared_ptr<PermutationGroup>(new PermutationGroup(
				ssc.construct(groupGenerators.begin(), groupGenerators.end(), groupBase.begin(), groupBase.end())
			));
		}
	} else {
		YALLOG_INFO(logger, "computation with empty symmetry group");
		pg = boost::shared_ptr<PermutationGroup>(new PermutationGroup(poly->realRowNumber()));
	}
	
	BOOST_ASSERT( pg );
    
	if (vm.count("automorphisms-only")) {
		if (yal::ReportLevel::get() >= yal::DEBUG3) {
			cout << "transversal sizes for Aut(P) " << endl;
			BOOST_FOREACH( const TRANSVERSAL &u, pg->U ) {
				cout << u.size() << " * ";
			}
			cout << endl;
		}
		cout << "generators for Aut(P) of order " << pg->order() << ":" << endl;
		BOOST_FOREACH( const PERMptr &p, pg->S ) {
			cout << " " << *p << endl;
		}
		cout << "SymPol format:" << endl;
		print_group(pg);
	} else if (vm.count("estimation-only")) {
		std::list<FaceWithData> dummy;
		const double est = rayComp->estimate(*poly, dummy);
		cout << "LRS estimation " << est << endl;
	} else {
		RecursionStrategy* rs = 0;
		bool withAdjacencies = vm.count("adjacencies");

		if (vm.count("direct")) {
			rs = new RecursionStrategyDirect();
			if (withAdjacencies)
				YALLOG_WARNING(logger, "Option --adjacencies requires ADM. Use one of the adm options instead of direct computation.")
		}	else if (vm.count("adm")) {
			rs = new RecursionStrategyADM(vm["adm"].as<double>());
		}	else if (vm.count("adm-dim")) {
			rs = new RecursionStrategyADMDimension(vm["adm-dim"].as<uint>());
		}	else if (vm.count("adm-incidence")) {
			rs = new RecursionStrategyADMIncidence(vm["adm-incidence"].as<uint>());
		}	else if (vm.count("idm-adm")) {
			const vector<double>& thresholds = vm["idm-adm"].as<vector<double> >();
			BOOST_ASSERT( thresholds.size() >= 2 );
			rs = new RecursionStrategyIDMADM(thresholds[0], thresholds[1]);
		} else if (vm.count("idm-adm-level")) {
			const vector<uint>& levels = vm["idm-adm-level"].as<vector<uint> >();
			BOOST_ASSERT( levels.size() >= 2 );
			rs = new RecursionStrategyIDMADMLevel(levels[0], levels[1]);
		} else if (vm.count("adm-idm-level")) {
			const vector<uint>& levels = vm["adm-idm-level"].as<vector<uint> >();
			BOOST_ASSERT( levels.size() >= 2 );
			rs = new RecursionStrategyADMIDMLevel(levels[1], levels[0]);
		}
		
		if (rs) {
			srand(time(NULL));
			FacesUpToSymmetryList rd(*pg, false, withAdjacencies);
			bool ret = false;
			
			if (vm.count("resume")) {
				rs->setDumpfile(vm["resume"].as<string>());
				ret = rs->resumeComputation(rayComp, *poly, *pg, rd);
			} else if (vm.count("dump")) {
				rs->setDumpfile(vm["dump"].as<string>());
				ret = rs->enumerateRaysUpToSymmetry(rayComp, *poly, *pg, rd);
			} else {
				ret = rs->enumerateRaysUpToSymmetry(rayComp, *poly, *pg, rd);
			}
			
			delete rs;
				
			if (ret) {
				if (poly->representation() == Polyhedron::H)
					cout << "V-representation" << endl;
				else if (poly->representation() == Polyhedron::V)
					cout << "H-representation" << endl;

				ulong dimension = poly->dimension();
				bool homogenized = false;
				ulong size = rd.size();
				// if polyhedron is V-representation, it is homogenized
				if (poly->representation() == Polyhedron::V) {
					homogenized = true;
					// homogenized dimension
					--dimension;
					// print skips cone apex
					--size;
				}
				cout << "* UP TO SYMMETRY" << endl << "begin" << endl << size << " " << dimension << " rational" << endl;
				PolyhedronIO::write(rd, homogenized, cout);
				cout << "end" << endl << "permutation group" << endl << "* order " << pg->order() << endl << "* w.r.t. to the original inequalities/vertices" << endl;
				print_group(pg);

				if (withAdjacencies) {
					ulong apexIndex = 0xffffffff;
					if (homogenized)
						apexIndex = static_cast<ulong>(rd.firstVertexIndex());

					cout << endl << "graph adjacencies {" << endl;
					BOOST_FOREACH(const FaceWithDataPtr& faceData, make_pair(rd.begin(), rd.end())) {
						BOOST_FOREACH(const FaceWithDataPtr& adjacent, faceData->adjacencies) {
							cout << "  " << correct_id(faceData->id, apexIndex)  << " -- " << correct_id(adjacent->id, apexIndex) << ";" << endl;
						}
					}
					cout << "}" << endl;
				}
			} else {
				YALLOG_ERROR(logger, "result of enumerateRaysUpToSymmetry invalid, aborting");
			}
		} else {
			YALLOG_ERROR(logger, "No recursion strategy selected, aborting");
		}
	}
    
	if (bTimeMeasurement) {
		timeVal = UsageStats::processTimeUser() - timeVal;
		cout << "elapsed time: " << timeVal << " seconds" << endl;
	}
	
	// clean up computation instruments
	rayComp->finish();
	delete rayComp;
	delete poly;
	
	PolyhedronDataStorage::cleanupStorage();
	
	return 0;
}
