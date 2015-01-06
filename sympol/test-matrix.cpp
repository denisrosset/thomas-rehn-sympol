#include "common.h"

//#include <predicate/matrix_automorphism_predicate.h>
//#include <search/partition/matrix_refinement1.h>
#include <search/partition/matrix_automorphism_search.h>
#include <permutation.h>
#include <transversal/explicit_transversal.h>
#include <construct/schreier_sims_construction.h>

#include <symmetric_group.h>

#include "zmatrix.h"
#include <iostream>

using namespace std;
using namespace permlib;
using namespace permlib::partition;

int main() {
// 	const ulong n = 5;
// 	const ulong k = 3;
// 	
// 	ZMatrix mat(n, k);
// 	mat.at(1,1) = 0;
// 	mat.at(2,2) = 0;
// 	mat.at(3,3) = 2;
// 	mat.at(4,4) = 2;
// 	
// 	mat.at(1,4) = 1;
// 	mat.at(4,1) = 1;
// 	mat.at(0,4) = 1;
// 	mat.at(4,0) = 1;
// 	mat.at(2,4) = 1;
// 	mat.at(4,2) = 1;
// 	mat.at(0,1) = 2;
// 	mat.at(1,0) = 2;
// 	mat.at(2,1) = 2;
// 	mat.at(1,2) = 2;
// 	mat.at(2,0) = 2;
// 	mat.at(0,2) = 2;

	const ulong n = 25;
	const ulong k = 5;
	srand(time(NULL));
	ZMatrix mat(n, k);
	for (uint i=0; i<n; ++i) {
		for (uint j=i; j<n; ++j) {
			mat.at(i,j) = randomInt(k);
			mat.at(j,i) = mat.at(i,j);
		}
	}
	
	for (uint i=0; i<n; ++i) {
		mat.at(i,0) = mat.at(i,1);
	}
	for (uint i=0; i<n; ++i) {
		mat.at(0,i) = mat.at(i,0);
	}
	mat.at(0,0) = mat.at(1,1);
	
	std::cout << "Matrix" << std::endl;
	mat.print();
	
	typedef Permutation PERM;
	typedef SymmetricGroup<PERM>::TRANS TRANS;
	typedef ExplicitTransversal<PERM> TRANSRET;
	
// 	PERMlist groupGenerators;
// 	boost::shared_ptr<PERM> gen1(new PERM(n, std::string("1 2 3")));
// 	groupGenerators.push_back(gen1);
// 	boost::shared_ptr<PERM> gen2(new PERM(n, std::string("3 4")));
// 	groupGenerators.push_back(gen2);
// 	boost::shared_ptr<PERM> gen3(new PERM(n, std::string("4 5")));
// 	groupGenerators.push_back(gen3);
	
	//SchreierSimsConstruction<PERM, TRANS> ssc(n);
	//BSGS<PERM,TRANS> bsgs(ssc.construct(groupGenerators.begin(), groupGenerators.end()));
	//MatrixAutomorphismSearch<PERM, TRANS> mas(bsgs, false);
	SymmetricGroup<PERM> s_n(n);
	MatrixAutomorphismSearch<PERM, TRANS, TRANSRET> mas(s_n, false);
	mas.construct(mat);
	
	BSGS<PERM,TRANSRET> K(n);
	mas.search(K);
	std::cout << K << std::endl;
	
	return 0;
}