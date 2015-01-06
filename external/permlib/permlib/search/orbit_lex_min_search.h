// ---------------------------------------------------------------------------
//
//  This file is part of PermLib.
//
// Copyright (c) 2009-2010 Thomas Rehn <thomas@carmen76.de>
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. The name of the author may not be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ---------------------------------------------------------------------------

#ifndef ORBIT_LEX_MIN_SEARCH_H_
#define ORBIT_LEX_MIN_SEARCH_H_

#include <permlib/predicate/pointwise_stabilizer_predicate.h>
#include <permlib/change/conjugating_base_change.h>
#include <permlib/change/random_base_transpose.h>
#include <permlib/search/dset.h>

#include <vector>
#include <limits>
#include <boost/dynamic_bitset.hpp>

namespace permlib {

/// algorithm to find the lexicographically minimal set in an orbit
/**
 * implements the algorithm ``Finding the Smallest Image of a Set''
 * by Steve Linton, 2004
 */
template<class BSGSIN>
class OrbitLexMinSearch {
public:
	/// constructor
	/**
	 * @param bsgs the group to build the orbit from
	 */
	OrbitLexMinSearch(const BSGSIN& bsgs)
		: m_bsgs(bsgs), m_cbc(bsgs), m_dsetAction(bsgs.n), m_orb(m_bsgs.n), m_orbVector(m_bsgs.n, 0), m_orbVectorIndex(0) {}

	/// searches the lexicographically minimal element of an orbit
	/**
	 * @param element one element of the orbit
	 * @param stabilizer setwise stabilizer of given set in the group; if present may speed up computations; parameter is currently ignored
	 * @return lexicographically smallest orbit element
	 */
	dset lexMin(const dset& element, const BSGSIN* stabilizer = NULL);

	/// compares two sets given as dsets lexicographically
	/**
	 * This is different to a < b as dynamic_bitsets because
	 * we are interested in the the representation of the bitsets
	 * as sets.
	 */
	static bool isLexSmaller(const dset& a, const dset& b);

private:
	BSGSIN m_bsgs;
	const BSGSIN* m_bsgsStabilizer;
	typedef typename BSGSIN::PERMtype PERM;
	typedef std::vector<boost::shared_ptr<PERM> > PERMvec;

	struct Candidate {
		dset D;
		dset J;

		Candidate(dset D) : D(D), J(D.size()) {}
		Candidate(dset D, dset J) : D(D), J(J) {}

		void print(const char* prefix) const {
			std::cout << prefix <<  ".J = " << J << "  ; " << prefix << ".D = " << D << std::endl;
		}
	};

	typedef Candidate* CandidatePtr;

	ConjugatingBaseChange<PERM, typename BSGSIN::TRANStype, RandomBaseTranspose<PERM, typename BSGSIN::TRANStype> > m_cbc;
	DSetAction<PERM> m_dsetAction;

	bool lexMin(uint i, uint k, const BSGSIN* stabilizer, const std::list<CandidatePtr>& candidates, std::list<CandidatePtr>& candidatesNext, dset& M_i, std::list<ulong>& base, PERMvec& S_i);
	/// finds the least element of an orbit of one number
	ulong orbMin(ulong element, const PERMvec& generators);

	/// given a set of elements, finds orbit representatives
	/**
	 * @param element the input set
	 * @param generators group generators for the orbit
	 * @return a set that contains one orbit representative for each orbit in $(element)
	 */
	dset* orbRepresentatives(dset element, const PERMlist& generators);

	// temporary variables for the orbMin calculation
	dset m_orb;
	std::vector<ulong> m_orbVector;
	uint m_orbVectorIndex;
};


template<class BSGSIN>
inline dset OrbitLexMinSearch<BSGSIN>::lexMin(const dset& element, const BSGSIN* stabilizer) {
	if (element.count() == element.size())
		return element;
	if (element.count() == 0)
		return element;
	CandidatePtr c0(new Candidate(element));

	std::list<CandidatePtr> candList0, candList1;
	std::list<CandidatePtr>* cand0 = &candList0;
	std::list<CandidatePtr>* cand1 = &candList1;

	cand0->push_back(c0);
	dset M_i(element.size());
	std::list<ulong> base;
	PERMvec S_i;
	S_i.reserve(m_bsgs.S.size());

	for (uint i = 0; i < element.count(); ++i) {
		if (lexMin(i, element.count(), stabilizer, *cand0, *cand1, M_i, base, S_i))
			break;
		std::swap(cand0, cand1);
	}
	std::for_each(candList0.begin(), candList0.end(), delete_object());
	std::for_each(candList1.begin(), candList1.end(), delete_object());

	return M_i;
}

template<class BSGSIN>
inline bool OrbitLexMinSearch<BSGSIN>::lexMin(uint i, uint k, const BSGSIN* stabilizer, const std::list<CandidatePtr>& candidates, std::list<CandidatePtr>& candidatesNext, dset& M_i, std::list<ulong>& base, PERMvec& S_i) {
	DEBUG(std::cout << "### START " << i << " with #" << candidates.size() << std::endl;)

	// if current stabilizer in the stabilizer chain is trivial we may
	// choose the minimal candidate and abort the search
	bool allOne = true;
	for (uint j = i; j < m_bsgs.B.size(); ++j) {
		if (m_bsgs.U[j].size() > 1) {
			allOne = false;
			break;
		}
	}
	if (allOne) {
		M_i = candidates.front()->D;
		BOOST_FOREACH(const CandidatePtr& R, candidates) {
			if (isLexSmaller(R->D, M_i)) {
				M_i = R->D;
			}
		}
		return true;
	}

	uint m = m_bsgs.n + 1;
	S_i.clear();
	PointwiseStabilizerPredicate<PERM> stab_i(m_bsgs.B.begin(), m_bsgs.B.begin() + i);
	std::copy_if(m_bsgs.S.begin(), m_bsgs.S.end(), std::back_inserter(S_i), stab_i);
	const ulong UNDEFINED_ORBIT = std::numeric_limits<ulong>::max();
	std::vector<ulong> orbitCache(m_bsgs.n, UNDEFINED_ORBIT);
	std::list<CandidatePtr> pass;

	BOOST_FOREACH(const CandidatePtr& R, candidates) {
		ulong m_R = m;
		for (ulong j = 0; j < R->D.size(); ++j) {
			if (R->J[j] || !R->D[j])
				continue;

			ulong val = orbitCache[j];
			if (val == UNDEFINED_ORBIT) {
				val = orbMin(j, S_i);
				orbitCache[j] = val;
			}
			if (m_R > val)
				m_R = val;
		}

		if (m_R < m) {
			m = m_R;
			pass.clear();
			pass.push_back(R);
		} else if (m_R == m) {
			pass.push_back(R);
		}
	}

	DEBUG(std::cout << " found m = " << m << std::endl;)
	M_i.set(m, 1);
	if (i == k-1)
		return true;

	base.push_back(m);
	m_cbc.change(m_bsgs, base.begin(), base.end());

	std::for_each(candidatesNext.begin(), candidatesNext.end(), delete_object());
	candidatesNext.clear();

	PERM* UNDEFINED_TRANSVERSAL = reinterpret_cast<PERM*>(1L);
	std::vector<PERM*> transversalCache(m_bsgs.n);
	BOOST_FOREACH(PERM*& p, transversalCache) {
		p = UNDEFINED_TRANSVERSAL;
	}
	BOOST_FOREACH(const CandidatePtr& R, pass) {
		for (ulong j = 0; j < R->D.size(); ++j) {
			if (!R->D[j])
				continue;

			PERM* perm = transversalCache[j];
			if (perm == UNDEFINED_TRANSVERSAL) {
				perm = m_bsgs.U[i].at(j);
				if (perm) {
					perm->invertInplace();
				}
				transversalCache[j] = perm;
			}

			if (!perm)
				continue;

			CandidatePtr c(new Candidate(R->D, R->J));
			m_dsetAction.apply(*perm, R->D, c->D);
			c->J.set(m);
			candidatesNext.push_back(c);
		}
	}

	BOOST_FOREACH(PERM* p, transversalCache) {
		if (p != UNDEFINED_TRANSVERSAL)
			delete p;
	}
	return false;
}

template<class BSGSIN>
inline ulong OrbitLexMinSearch<BSGSIN>::orbMin(ulong element, const PERMvec& generators) {
	if (element == 0)
		return 0;

	ulong minElement = element;
	m_orb.reset();
	m_orb.set(element, 1);
	m_orbVectorIndex = 0;
	m_orbVector[m_orbVectorIndex++] = element;

	for (uint i = 0; i < m_orbVectorIndex; ++i) {
		const ulong &alpha = m_orbVector[i];
		BOOST_FOREACH(const PERMptr& p, generators) {
			ulong alpha_p = *p / alpha;
			if (alpha_p == 0)
				return 0;
			if (!m_orb[alpha_p]) {
				m_orbVector[m_orbVectorIndex++] = alpha_p;
				m_orb.set(alpha_p);
				if (alpha_p < minElement)
					minElement = alpha_p;
			}
		}
	}

	return minElement;
}


template<class BSGSIN>
inline dset* OrbitLexMinSearch<BSGSIN>::orbRepresentatives(dset element, const PERMlist& generators) {
	dset* ret = new dset(element.size());

	for (uint j = 0; j < element.size(); ++j) {
		if (!element[j])
			continue;

		m_orb.set();
		m_orb.set(j, 0);
		m_orbVectorIndex = 0;
		m_orbVector[m_orbVectorIndex++] = j;
		for (uint i = 0; i < m_orbVectorIndex; ++i) {
			const ulong &alpha = m_orbVector[i];
			BOOST_FOREACH(const PERMptr& p, generators) {
				ulong alpha_p = *p / alpha;
				if (m_orb[alpha_p]) {
					m_orbVector[m_orbVectorIndex++] = alpha_p;
					m_orb.reset(alpha_p);
				}
			}
		}

		element &= m_orb;
		ret->set(j);
	}

	return ret;
}


template<class BSGSIN>
inline bool OrbitLexMinSearch<BSGSIN>::isLexSmaller(const dset& a, const dset& b) {
		dset::size_type i = a.find_first(), j = b.find_first();
		while (i != dset::npos && j != dset::npos) {
			if (i < j)
				return true;
			if (i > j)
				return false;
			i = a.find_next(i);
			j = b.find_next(j);
		}
		return false;
	}

} // ns permlib

#endif // ORBIT_LEX_MIN_SEARCH_H_
