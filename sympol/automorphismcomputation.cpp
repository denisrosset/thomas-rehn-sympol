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

#include "automorphismcomputation.h"

#if HAVE_NAUTY && HAVE_NTL
extern "C" {
  #include <nauty.h>
}
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#endif // NAUTY && NTL

#include <vector>
#include <permlib/construct/schreier_sims_construction.h>
#include <permlib/search/partition/matrix_automorphism_search.h>
#include <permlib/symmetric_group.h>

#include "matrix/zmatrix.h"
#include "matrix/invert.h"

using namespace sympol;
using namespace sympol::matrix;
using namespace permlib;
using namespace yal;

LoggerPtr AutomorphismComputation::logger(Logger::getLogger("AutomComp "));

boost::shared_ptr< PermutationGroup > AutomorphismComputation::computeRestrictedIsomorphisms ( const Polyhedron& poly ) {
	// we assume that the polyhedron is a full-dimensional cone
	// so the rank of the inequality matrix is one less than the 
	// inequality vectors have entries;
	const ulong matrixRank = poly.dimension() - 1;
	const ulong matrixRows = poly.rows();
    
	typedef std::map<mpq_class,uint> WeightMap;
	WeightMap weights;
	mpq_t temp;
	mpq_init(temp);
	
	//
	// calculate vertex weights
	//
	typedef Matrix<mpq_class> QMatrix;
	QMatrix Q(matrixRank);
	for (ulong i=0; i<matrixRank; ++i) {
		for (ulong j=0; j<matrixRank; ++j) {
			BOOST_FOREACH(const QArray& row, poly.rowPair()) {
				mpq_mul(temp, row[i+1], row[j+1]);
				mpq_add(Q.at(i,j).get_mpq_t(), Q.at(i,j).get_mpq_t(), temp);
				//Q.at(i,j) += row[i] * row[j];
			}
		}
	}
	YALLOG_DEBUG2(logger, "Q = " << std::endl << Q);
	
	QMatrix Qinv(matrixRank);
	if (!Invert<QMatrix>(&Q).invert(&Qinv)) {
		YALLOG_ERROR(logger, "could not invert matrix");
		return boost::shared_ptr< PermutationGroup >();
	}
	
	YALLOG_INFO(logger, "matrix inversion complete");
	YALLOG_DEBUG3(logger, "Qinv = " << std::endl << Qinv);
    
	ZMatrix zMatrix(matrixRows);
    
	uint weightIndex = 0;
	ulong i = 0, j = 0;
    BOOST_FOREACH(const QArray& row1, poly.rowPair()) {
		j = 0;
		BOOST_FOREACH(const QArray& row2, poly.rowPair()) {
			if (i < j)
				break;
			mpq_class newWeight;
			for (uint k = 0; k < matrixRank; ++k) {
				for (uint l = 0; l < matrixRank; ++l) {
					mpq_mul(temp, row1[k+1], row2[l+1]);
					mpq_mul(temp, temp, Qinv.at(k,l).get_mpq_t());
					mpq_add(newWeight.get_mpq_t(), newWeight.get_mpq_t(), temp);
					//newWeight += Qinv.at(k,l) * row1[k] * row2[l];
				}
			}
			std::pair<WeightMap::iterator, bool> suc = weights.insert(std::make_pair(newWeight, weightIndex));
			zMatrix.at(i,j) = suc.first->second;
			zMatrix.at(j,i) = suc.first->second;
			if (suc.second)
				++weightIndex;
			++j;
		}
		++i;
	}
	mpq_clear(temp);
	
	zMatrix.k() = weightIndex;
	
	YALLOG_DEBUG(logger, "zMatrix = " << std::endl << zMatrix);
	
	SymmetricGroup<PERM> s_n(matrixRows);
	partition::MatrixAutomorphismSearch<SymmetricGroup<PERM>, TRANSVERSAL> mas(s_n, false);
	mas.construct(zMatrix);
	
	BSGS<PERM,TRANSVERSAL>* K = new BSGS<PERM,TRANSVERSAL>(matrixRows);
	mas.search(*K);
	YALLOG_INFO(logger, "matrix automorphism search complete; found group of order " << K->order());
	
	return boost::shared_ptr<PermutationGroup>(K);
}


#if HAVE_NAUTY && HAVE_NTL
// some global variables for use with nauty callback nauty
static std::list<boost::shared_ptr<Permutation> > s_generators;
static ulong permN = 0L;

// nauty callback function
static void setAutomorphisms(int count, permutation* perm, int* orbits, int numorbits, int stabvertex, int n);

boost::shared_ptr< PermutationGroup > AutomorphismComputation::computeRestrictedIsomorphismsNauty ( const Polyhedron& poly ) {
    // we assume that the polyhedron is a full-dimensional cone
    // so the rank of the inequality matrix is one less than the 
    // inequality vectors have entries;
    const ulong matrixRank = poly.dimension() - 1;
    const ulong matrixRows = poly.rows();
    
    //
    // calculate vertex weights
    //
    NTL::mat_ZZ Q;
    NTL::vec_ZZ* v = new NTL::vec_ZZ[matrixRows];
    NTL::mat_ZZ Qinv;
    Q.SetDims(matrixRank, matrixRank);
    NTL::ZZ acc, tmp, c1, c2;
    char *s;
    std::set<NTL::ZZ> weights;
    std::vector<NTL::ZZ> vecWeights;
    const NTL::ZZ** zWeights = new const NTL::ZZ*[matrixRows*matrixRows];

    // for converting rational coefficients into integer
    mpz_t lcm;
    mpz_t gmp_tmp;
    mpz_init(lcm);
    mpz_init(gmp_tmp);

    for (ulong i=0; i<matrixRows; ++i) {
        v[i].SetLength(matrixRank);
    }
    ulong jIndex = 0;
    for (ulong j=0; j<matrixRank; ++j) {
        ulong i = 0;  
        for (Polyhedron::RowIterator it = poly.rowsBegin(); it != poly.rowsEnd(); ++it) {
            // NTL wants integer coefficients, so we have to convert them, if we have some rational ones
            // TODO: for performance reasons we could cache the lcm value per row/i
            (*it).denominatorLCM(lcm);
            mpz_mul(gmp_tmp, lcm, mpq_numref(poly.element(i, j+1)));
            // TODO: make this code nicer
            // hack to convert GMP numbers into NTL numbers
            // number system is decimal (10)
            s = mpz_get_str (NULL, 10, gmp_tmp);
            v[i][jIndex] = NTL::to_ZZ (s);
            free (s);
            ++i;
        }
        ++jIndex;
    }
    
    // we don't need the variable any more
    mpz_clear(lcm);
    mpz_clear(gmp_tmp);
    
    for (long i=0; i<Q.NumRows(); ++i) {
        for (long j=0; j<Q.NumCols(); ++j) {
            acc = 0L;
            for (ulong m=0; m<matrixRows; ++m) {
                NTL::mul(tmp, v[m][i], v[m][j]);
                NTL::add(acc, acc, tmp);
            }
            Q[i][j] = acc;
        }
    }
    
    NTL::inv(tmp, Qinv, Q, 1);
    
    NTL::vec_ZZ r1;
    for (ulong i=0; i<matrixRows; ++i) {
        for (ulong j=i; j<matrixRows; ++j) {
            NTL::mul(r1, Qinv, v[j]);
            NTL::InnerProduct(tmp, v[i], r1);

            NTL::ZZ newWeight(tmp);
            std::pair<std::set<NTL::ZZ>::iterator, bool> suc = weights.insert(newWeight);
            const NTL::ZZ* pNewWeight = &(*(suc.first));
            zWeights[i*matrixRows + j] = pNewWeight;
            zWeights[j*matrixRows + i] = pNewWeight;
        }
    }

    std::set<NTL::ZZ>::const_iterator it;
    int r = 0;
    for (it = weights.begin(); it != weights.end(); ++it) {
        vecWeights.push_back(*it);
        ++r;
    }
    
    const ulong ulK = weights.size();
    //
    // instrumentalize nauty
    //
    
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk(stats);

    const ulong ulKBits = (ulong)(log2(ulK)+1.0);
    int n,m;
    set *gv;
    n = matrixRows * ulKBits;
    //options.writeautoms = TRUE;
    options.defaultptn = FALSE;
    // WARNING: NOT THREAD SAFE - static global variable for nauty-C
    s_generators.clear();
    permN = matrixRows;
    // set nauty callback function
    options.userautomproc = &setAutomorphisms;
    
    m = (n + WORDSIZE - 1) / WORDSIZE;
    nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
    // 50*m recommended workspace setting in nauty manual
    const int worksize = 50 * m;
    setword *workspace = new setword[worksize];
    
    char* sMalloc = "malloc";
    DYNALLOC2(graph,g,g_sz,m,n,sMalloc);
    DYNALLOC1(int,lab,lab_sz,n,sMalloc);
    DYNALLOC1(int,ptn,ptn_sz,n,sMalloc);
    DYNALLOC1(int,orbits,orbits_sz,n,sMalloc);

    // convert from edge-weighted to vertex-weighted graph
    // for automorphism group calculation
    for (int _v = 0; _v < n; ++_v) {
        gv = GRAPHROW(g,_v,m);

        const ulong kIndex = _v % ulKBits;
        const ulong nIndex = _v / ulKBits;
        EMPTYSET(gv,m);
        for (ulong w = nIndex*ulKBits; w < (nIndex+1)*ulKBits; ++w) {
            ADDELEMENT(gv,w);
        }
        
        for (ulong w=0; w<matrixRows; ++w) {
            // determine color index kb
            for (ulong kb=0; kb < ulK; ++kb) {
                if (vecWeights[kb] == *(zWeights[nIndex*matrixRows + w])) {
                    // check, if bit is set and an edge has to be inserted
                    // cf. nauty user manual 2.4b3, p.25
                    if ((1 << kIndex) & kb) {
                        ADDELEMENT(gv,w * ulKBits + kIndex);
                    }
                }
            }
        }
    }
    ulong _c = 0;
    for (ulong c=0; c<ulKBits; ++c) {
        for (ulong _v = 0; _v < matrixRows; ++_v) {
            lab[_c] = _v * ulKBits + c;
            ptn[_c] = NAUTY_INFINITY;
            ++_c;
        }
        ptn[_c-1] = 0;
    }

    // call nauty for automorphism calculation
    nauty(g, lab, ptn, NULL, orbits, &options, &stats,
        workspace, worksize, m, n, NULL);
    
    delete[] workspace;
    delete[] v;
    delete[] zWeights;
  
  SchreierSimsConstruction<PERM, TRANSVERSAL> schreierSims(permN);
  return boost::shared_ptr<PermutationGroup>(new PermutationGroup(schreierSims.construct(s_generators.begin(),
                              s_generators.end())));
}

/**
 * nauty callback function according to nauty API
 * called on each automorphism found by nauty
 */
static void setAutomorphisms(int count, permutation* perm, int* orbits, int numorbits, int stabvertex, int n) {
    // transform vertex-weighted graph back into edge-weighted
    std::vector<ulong> filteredPerm(permN);
    // k is the number in each vertex-equivalence class
    int k = n / permN;
    for (int i=0; i<n; i+=k) {
        filteredPerm[i/k] = perm[i] / k;
    }
    boost::shared_ptr<Permutation> gen(new Permutation(filteredPerm));
  s_generators.push_back(gen);
    //std::cout << "Perm: " << *gen << std::endl;
  //permlib::print_iterable(filteredPerm.begin(), filteredPerm.end(), 0, "PERM");
}
#endif
