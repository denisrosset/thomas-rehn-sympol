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


#ifndef BASE_SEARCH_H_
#define BASE_SEARCH_H_

#include <permlib/change/conjugating_base_change.h>
#include <permlib/change/random_base_transpose.h>

#include <boost/dynamic_bitset.hpp>

namespace permlib {

/// base class for searching in a group
template<class BSGSIN, class TRANSRET>
class BaseSearch {
public:
	typedef typename BSGSIN::PERMtype PERM;
	typedef typename BSGSIN::TRANStype TRANS;
	typedef boost::shared_ptr< typename BaseSearch<BSGSIN,TRANSRET>::PERM> PERMptrType;
	typedef std::list<PERMptrType> PERMlistType;
	
	/// constructor
	/**
	 * @param bsgs BSGS to search in
	 * @param pruningLevelDCM prune levels smaller than pruningLevelDCM by double coset minimality with base change
	 * @param stopAfterFirstElement true iff the search can be stopped after the first element found with the desired property
	 */
	BaseSearch(const BSGSIN& bsgs, uint pruningLevelDCM, bool stopAfterFirstElement);
	
	/// destructor
	virtual ~BaseSearch(){}
	
	/// finds minimal elements in an orbit
	/**
	 * returns true iff beta_i is the minimal element of the orbit of alpha under action of the i-th stabilizer of groupK
	 */
	bool minOrbit(ulong alpha, BSGS<PERM,TRANSRET> &groupK, uint i, ulong beta_i) const;
	
	/// searches for a coset representative if one exists
	virtual PERMptrType searchCosetRepresentative();
	
	/// searches for a coset representative if one exists
	/**
	 * the two arguments are two groups K and L such that \f$KgL \subseteq G(\mathcal P)\quad \Leftrightarrow g \in G(\mathcal P)\f$
	 * @param groupK subgroup of G
	 * @param groupL subgroup of G
	 * @return pointer to a coset representative element or NULL
	 */
	virtual PERMptrType searchCosetRepresentative(BSGS<PERM,TRANSRET> &groupK, BSGS<PERM,TRANSRET> &groupL) = 0;
	
	/// nodes visited during backtrack search
	ulong m_statNodesVisited;
	/// number of nodes where (simple) double coset minimality pruning was in effect
	ulong m_statNodesPrunedCosetMinimality;
	/// number of nodes where advanced double coset minimality pruning with base change was in effect
	ulong m_statNodesPrunedCosetMinimality2;
	/// number of nodes where a child constraint pruning was in effect
	ulong m_statNodesPrunedChildRestriction;

protected:
	/// main BSGS to search in
	BSGSIN m_bsgs;
	/// second BSGS of a group the sough elements have to member of
	BSGSIN* m_bsgs2;
	/// predicate that matches sought elements
	boost::scoped_ptr<SubgroupPredicate<PERM> > m_pred;
	
	/// base point order
	std::vector<ulong> m_order;
	/// a sorter with respect to m_order
	boost::scoped_ptr<BaseSorterByReference> m_sorter;
	
	/// base change algorithm
	ConjugatingBaseChange<PERM,TRANS,RandomBaseTranspose<PERM,TRANS> > m_baseChange;

	/// leves i with 0 <= i < m_pruningLevelDCM are prunged by advanced double coset minimality 
	const uint m_pruningLevelDCM;
	/// true iff other m_limit variables have been initialized
	bool m_limitInitialized;
	/// number of base points that correspond to maximal backtrack level m_limitLevel
	uint m_limitBase;
	/// maximal backtrack level
	uint m_limitLevel;
	
	/// try to prune with advanced double coset minimality
	bool pruneDCM(const PERM& t, uint backtrackLevel, BSGS<PERM,TRANSRET> &groupK, BSGS<PERM,TRANSRET> &groupL);
	/// true iff level is a leaf level
	bool checkLeaf(uint level);
	/// processes a leaf and adds corresponding element to the generator set of K
	uint processLeaf(const PERM& t, uint level, uint backtrackLevel, uint completed, BSGS<PERM,TRANSRET> &groupK, BSGS<PERM,TRANSRET> &groupL);
	/// base of the sought subgroup
	virtual const std::vector<ulong>& subgroupBase() const = 0;
	
	/// sets up a BSGS structure for an empty group with base subgroupBase()
	void setupEmptySubgroup(BSGS<PERM,TRANSRET>& group) const;
	
	/// true iff the search can be stopped after the first element found with the desired property
	const bool m_stopAfterFirstElement;
	/// last element found with desired property; only used if m_stopAfterFirstElement is true
	PERMptrType m_lastElement;
private:
	static PERMlistType ms_emptyList;
};

//
// IMPLEMENTATION
//

template<class BSGSIN,class TRANSRET>
typename BaseSearch<BSGSIN,TRANSRET>::PERMlistType BaseSearch<BSGSIN,TRANSRET>::ms_emptyList;


template<class BSGSIN,class TRANSRET>
BaseSearch<BSGSIN,TRANSRET>::BaseSearch(const BSGSIN& bsgs, uint pruningLevelDCM, bool stopAfterFirstElement) 
	: m_statNodesVisited(0), m_statNodesPrunedCosetMinimality(0), m_statNodesPrunedCosetMinimality2(0),
	  m_statNodesPrunedChildRestriction(0),
	  m_bsgs(bsgs), m_bsgs2(0), m_pred(0), m_baseChange(m_bsgs),
	  m_pruningLevelDCM(pruningLevelDCM),
	  m_limitInitialized(false), m_limitBase(0), m_limitLevel(0),
	  m_stopAfterFirstElement(stopAfterFirstElement),
	  m_lastElement()
{ 
}


template<class BSGSIN,class TRANSRET>
bool BaseSearch<BSGSIN,TRANSRET>::minOrbit(ulong alpha, BSGS<PERM,TRANSRET> &groupK, uint i, ulong beta_i) const {
	PERMlistType S_i;
	std::copy_if(groupK.S.begin(), groupK.S.end(), std::back_inserter(S_i), PointwiseStabilizerPredicate<PERM>(groupK.B.begin(), groupK.B.begin() + i));
	if (S_i.empty()) {
		if (alpha == beta_i)
			return true;
		return (*m_sorter)(beta_i, alpha);
	}
	
	//TODO: avoid multiple allocation?
	boost::dynamic_bitset<> orbitCharacteristic(m_bsgs.n);
	orbitCharacteristic.set(alpha, 1);
	std::list<ulong> orbit;
	orbit.push_back(alpha);
    for (std::list<ulong>::const_iterator it = orbit.begin(); it != orbit.end(); ++it) {
        ulong beta = *it;
        BOOST_FOREACH(const PERMptrType& p, S_i) {
            ulong beta_p = *p / beta;
            if (!orbitCharacteristic[beta_p]) {
                orbitCharacteristic.set(beta_p, 1);
                orbit.push_back(beta_p);
				if ((*m_sorter)(beta_p, beta_i)) {
					DEBUG(std::cout << "DCM2 beta_p = " << beta_p+1 << " , beta_i = " << beta_i+1 << std::endl;)
					return false;
				}
            }
        }
    }
    return true;
}

template<class BSGSIN,class TRANSRET>
bool BaseSearch<BSGSIN,TRANSRET>::pruneDCM(const PERM& t, uint backtrackLevel, BSGS<PERM,TRANSRET> &groupK, BSGS<PERM,TRANSRET> &groupL) {
	// change base only for the lower nodes in the tree
	if (backtrackLevel < m_pruningLevelDCM) {
		//TODO: avoid multiple allocation?
		std::vector<ulong> newBaseImage(subgroupBase().begin(), subgroupBase().end());
		for (uint j=0; j<=backtrackLevel; ++j)
			newBaseImage[j] = t / newBaseImage[j];
		//print_iterable(newBaseImage.begin(), newBaseImage.begin() + (backtrackLevel+1), 1, "new base image");
		ConjugatingBaseChange<PERM,TRANSRET,RandomBaseTranspose<PERM,TRANSRET> > cbc(groupL);
		cbc.change(groupL, newBaseImage.begin(), newBaseImage.begin() + (backtrackLevel+1), false);
		//print_iterable(groupL.B.begin(), groupL.B.end(), 1, "new base");
	}

	const ulong alpha = groupK.B[backtrackLevel];
	
	for (uint i = 0; i <= backtrackLevel; ++i) {
		if (i == backtrackLevel || groupK.U[i].contains(alpha)) {
			DEBUG(std::cout << "DCM2 found " << (alpha+1) << " in U_" << i << "  btLevel " << backtrackLevel << std::endl;)
			DEBUG(std::cout << "     t = " << t << std::endl;)
			
			if (!minOrbit(t / alpha, groupL, i, t / groupK.B[i])) {
				DEBUG(std::cout << "DCM2 : " << ((t / groupK.B[i]) + 1) << " // " << ((t / alpha) + 1) << std::endl;)
				DEBUG(std::cout << " K = " << groupK << std::endl;)
				DEBUG(std::cout << " L = " << groupL << std::endl;)
				return true;
			}
		}
		if (t / groupK.B[i] != groupL.B[i])
			return false;
	}
	return false;
}

template<class BSGSIN,class TRANSRET>
bool BaseSearch<BSGSIN,TRANSRET>::checkLeaf(uint level) {
	return m_limitInitialized && level >= m_limitLevel;
}

template<class BSGSIN,class TRANSRET>
uint BaseSearch<BSGSIN,TRANSRET>::processLeaf(const PERM& t, uint level, uint backtrackLevel, uint completed, BSGS<PERM,TRANSRET> &groupK, BSGS<PERM,TRANSRET> &groupL) {
	DEBUG(std::cout << "XXX level " << level << "  bLevel " << backtrackLevel << std::endl;)
	DEBUG(std::cout << "XXX limitLevel " << m_limitLevel << "  limitBase " << m_limitBase << std::endl;)
	if ((*m_pred)(t)) {
		if (m_stopAfterFirstElement) {
			m_lastElement = PERMptrType(new PERM(t));
			return 0;
		}
		const bool isIdentity = t.isIdentity();
		int posK = 0, posL = 0;
		if (m_limitInitialized && level == m_limitLevel && isIdentity) {
			PointwiseStabilizerPredicate<PERM> stabPred(m_bsgs.B.begin(), m_bsgs.B.begin() + m_limitBase);
			BOOST_FOREACH(const PERMptrType &s, m_bsgs.S) {
				if (stabPred(s)) {
					DEBUG(std::cout << *s << " extended gen\n";)
					BOOST_ASSERT((*m_pred)(*s));
					PERMptrType sK(new PERM(*s));
					PERMptrType sL(new PERM(*s));
					posK = std::max(posK, groupK.insertGenerator(sK, false));
					posL = std::max(posL, groupL.insertGenerator(sL, false));
				}
			}
			//return completed;
		}
		if (!isIdentity) {
			PERMptrType genK(new PERM(t));
			posK = std::max(posK, groupK.insertGenerator(genK, false));
			PERMptrType genL(new PERM(t));
			posL = std::max(posL, groupL.insertGenerator(genL, false));
			DEBUG(std::cout << "-- accepted" << std::endl;)
		}
		groupK.updateOrbits(posK);
		groupL.updateOrbits(posL);
		return completed;
	}
	return level;
}

template<class BSGSIN,class TRANSRET>
void BaseSearch<BSGSIN,TRANSRET>:: setupEmptySubgroup(BSGS<PERM,TRANSRET>& group) const {
	group.B = subgroupBase();
	group.U.resize(subgroupBase().size(), TRANSRET(this->m_bsgs.n));
	for (uint i=0; i<subgroupBase().size(); ++i)
		group.orbit(i, ms_emptyList);
}

template<class BSGSIN,class TRANSRET>
typename BaseSearch<BSGSIN,TRANSRET>::PERMptrType BaseSearch<BSGSIN,TRANSRET>::searchCosetRepresentative() {
	BSGS<PERM,TRANSRET> groupK(this->m_bsgs.n);
	BSGS<PERM,TRANSRET> groupL(this->m_bsgs.n);
	
	setupEmptySubgroup(groupK);
	setupEmptySubgroup(groupL);
	
	return this->searchCosetRepresentative(groupK, groupL);
}

}

#endif // -- BASE_SEARCH_H_
