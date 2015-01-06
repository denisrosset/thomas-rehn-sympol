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


#ifndef RBASE_H_
#define RBASE_H_

#include <permlib/predicate/subgroup_predicate.h>

#include <permlib/search/base_search.h>

#include <permlib/search/partition/partition.h>
#include <permlib/search/partition/refinement_family.h>
#include <permlib/search/partition/backtrack_refinement.h>

#include <permlib/change/conjugating_base_change.h>
#include <permlib/change/random_base_transpose.h>

#include <permlib/sorter/base_sorter.h>

#include <utility>

namespace permlib {
namespace partition {

/// R-base for partition backtracking
template<class BSGSIN,class TRANSRET>
class RBase : public BaseSearch<BSGSIN,TRANSRET> {
public:
	typedef typename BaseSearch<BSGSIN,TRANSRET>::PERM PERM;
	typedef typename BaseSearch<BSGSIN,TRANSRET>::TRANS TRANS;
	typedef boost::shared_ptr< typename RBase<BSGSIN,TRANSRET>::PERM> PERMptrType;
	
	/// constructor
	/**
	 * @param bsgs BSGS to search in
	 * @param pruningLevelDCM prune levels smaller than pruningLevelDCM by double coset minimality with base change
	 * @param stopAfterFirstElement true iff the search can be stopped after the first element found with the desired property
	 */
	RBase(const BSGSIN& bsgs, uint pruningLevelDCM, bool stopAfterFirstElement = false);
	
	typedef typename Refinement<PERM>::RefinementPtr RefinementPtr;
	typedef typename RefinementFamily<PERM>::PartitionPtr PartitionPtr;
	typedef typename std::list<std::pair<PartitionPtr,RefinementPtr> >::const_iterator PartitionIt;
	
	/// perform search and store result in groupK
	void search(BSGS<PERM,TRANSRET> &groupK);
	
	using BaseSearch<BSGSIN,TRANSRET>::searchCosetRepresentative;
	virtual PERMptrType searchCosetRepresentative(BSGS<PERM,TRANSRET> &groupK, BSGS<PERM,TRANSRET> &groupL);
protected:
	/// partition to base the backtrack tree on
	Partition m_partition;
	Partition m_partition2;
	
	/// constructs an R-base for given predicate and refinement family
	/**
	 * group membership \f$\mathcal P\f$-refinement for m_bsgs is used by default
	 * @param pred 
	 * @param predRefinement refinement family to use; may be zero to use only group membership \f$\mathcal P\f$-refinement
	 */
	void construct(SubgroupPredicate<PERM>* pred, RefinementFamily<PERM>* predRefinement);
	
	/// callback when a new fix point appears during R-base construction
	virtual uint processNewFixPoints(const Partition& pi, uint level);
	
	virtual const std::vector<ulong>& subgroupBase() const;
private:
	/// base of the sough subgroup based on R-base
	std::vector<ulong> m_subgroupBase;
	/// actual R-base
	std::list<std::pair<PartitionPtr,RefinementPtr> > partitions;
	
	/// recursive backtrack search
	uint search(PartitionIt pIt, Partition &pi, const PERM& t, const PERM* t2, uint level, uint backtrackLevel, uint& completed, BSGS<PERM,TRANSRET> &groupK, BSGS<PERM,TRANSRET> &groupL);
	
	/// updates t2 such that sigma^t2 = pi2
	bool updateMappingPermutation(const BSGSIN& bsgs, const Partition& sigma, const Partition& pi2, PERM& t2) const;
};

template<class BSGSIN,class TRANSRET>
RBase<BSGSIN,TRANSRET>::RBase(const BSGSIN& bsgs, uint pruningLevelDCM, bool stopAfterFirstElement) 
	: BaseSearch<BSGSIN,TRANSRET>(bsgs, pruningLevelDCM, stopAfterFirstElement),
	  m_partition(bsgs.n), m_partition2(bsgs.n)
{ }

template<class BSGSIN,class TRANSRET>
void RBase<BSGSIN,TRANSRET>::construct(SubgroupPredicate<PERM>* pred, RefinementFamily<PERM>* predRefinement) {
	this->m_pred.reset(pred);
	typedef typename boost::shared_ptr<RefinementFamily<PERM> > RefinementFamilyPtr;
	std::list<RefinementFamilyPtr> refinements;
	
	if (!this->m_bsgs.isSymmetricGroup()) {
		RefinementFamilyPtr gr( new GroupRefinementFamily<PERM,TRANS>(this->m_bsgs) );
		refinements.push_back(gr);
	}
	
	if (predRefinement) {
		RefinementFamilyPtr predR( predRefinement );
		refinements.push_back(predR);
	}
	
	DEBUG(print_iterable(this->m_bsgs.B.begin(), this->m_bsgs.B.end(), 1, "orig BSGS");)
	
	Partition pi(m_partition);
	while (pi.cells() < this->m_bsgs.n) {
		DEBUG(std::cout << std::endl << "PI0 = " << pi << std::endl;)
		bool found = false;
		do {
			found = false;
			uint foo = 0;
			BOOST_FOREACH(RefinementFamilyPtr ref, refinements) {
				const uint oldFixPointsSize = pi.fixPointsSize();
				std::pair<PartitionPtr,RefinementPtr> newRef = ref->apply(pi);
				if (newRef.first) {
					partitions.push_back(newRef);
					if (oldFixPointsSize < pi.fixPointsSize()) {
						processNewFixPoints(pi, partitions.size());
					}
					//std::cout << "BSGS " << this->m_bsgs;
					found = true;
				}
				++foo;
			}
		} while(found);
		
		DEBUG(std::cout << std::endl << "PI1 = " << pi << std::endl;)
		
		if (pi.cells() < this->m_bsgs.n) {
			ulong alpha = -1;
			//print_iterable(pi.fixPointsBegin(), pi.fixPointsEnd(), 1, "  fix0");
			//print_iterable(this->m_bsgs.B.begin(), this->m_bsgs.B.end(), 1, "bsgs0");
			if (pi.fixPointsSize() < this->m_bsgs.B.size())
				alpha = this->m_bsgs.B[pi.fixPointsSize()];
			DEBUG(std::cout << "choose alpha = " << alpha << std::endl;)
			RefinementPtr br(new BacktrackRefinement<PERM>(this->m_bsgs.n, alpha));
			BacktrackRefinement<PERM>* ref = dynamic_cast<BacktrackRefinement<PERM> *>(br.get());
			ref->initializeAndApply(pi);
			PartitionPtr newPi(new Partition(pi));
			DEBUG(std::cout << "BACKTRACK " << (ref->alpha()+1) << " in " << pi << "    -->    " << *newPi << std::endl;)
			partitions.push_back(std::make_pair(newPi, br));
			
			processNewFixPoints(pi, partitions.size());
			
			//std::cout << "BSGS " << this->m_bsgs;
			m_subgroupBase.push_back(ref->alpha());
		}
	}
	
	this->m_order = BaseSorterByReference::createOrder(this->m_bsgs.n, pi.fixPointsBegin(), pi.fixPointsEnd());
	this->m_sorter.reset(new BaseSorterByReference(this->m_order));
	for (typename std::list<std::pair<PartitionPtr,RefinementPtr> >::iterator pIt = partitions.begin(); pIt != partitions.end(); ++pIt) {
		(*pIt).second->sort(*this->m_sorter, 0);
		DEBUG(std::cout << "SIGMA = " << *(*pIt).first << std::endl;)
	}
	
	DEBUG(print_iterable(this->m_order.begin(), this->m_order.end(), 0, "ORDER");)
}

template<class BSGSIN,class TRANSRET>
uint RBase<BSGSIN,TRANSRET>::processNewFixPoints(const Partition& pi, uint level) {
	const uint basePos = this->m_baseChange.change(this->m_bsgs, pi.fixPointsBegin(), pi.fixPointsEnd(), true);
	if (this->m_bsgs2)
		this->m_baseChange.change(*this->m_bsgs2, pi.fixPointsBegin(), pi.fixPointsEnd(), true);
	//print_iterable(pi.fixPointsBegin(), pi.fixPointsEnd(), 1, "  fix");
	DEBUG(print_iterable(this->m_bsgs.B.begin(), this->m_bsgs.B.end(), 1, "change base");)
	return basePos;
}

template<class BSGSIN,class TRANSRET>
void RBase<BSGSIN,TRANSRET>::search(BSGS<PERM,TRANSRET> &groupK) {
	BOOST_ASSERT( this->m_pred != 0 );
	
	setupEmptySubgroup(groupK);
	
	uint completed = partitions.size();
	BSGS<PERM,TRANSRET> groupL(groupK);
	PERM identH(this->m_bsgs.n);
	search(partitions.begin(), m_partition2, PERM(this->m_bsgs.n), &identH, 0, 0, completed, groupK, groupL);
}

template<class BSGSIN,class TRANSRET>
typename RBase<BSGSIN,TRANSRET>::PERMptrType RBase<BSGSIN,TRANSRET>::searchCosetRepresentative(BSGS<PERM,TRANSRET> &groupK, BSGS<PERM,TRANSRET> &groupL) {
	BOOST_ASSERT( this->m_pred != 0 );
	
	// !!!
	//
	//  TODO:  check that groupK and groupL have the right base (starting with subgroupBase)
	//
	// !!!
	
	uint completed = partitions.size();
	//BSGS<PERM,TRANS> groupL(groupK);
	PERM identH(this->m_bsgs.n);
	search(partitions.begin(), m_partition2, PERM(this->m_bsgs.n), &identH, 0, 0, completed, groupK, groupL);
	
	return BaseSearch<BSGSIN,TRANSRET>::m_lastElement;
}



template<class BSGSIN,class TRANSRET>
uint RBase<BSGSIN,TRANSRET>::search(PartitionIt pIt, Partition &pi, const PERM& t, const PERM* t2, uint level, uint backtrackLevel, uint& completed, BSGS<PERM,TRANSRET> &groupK, BSGS<PERM,TRANSRET> &groupL) {
	++this->m_statNodesVisited;

	if (pIt == partitions.end() || this->checkLeaf(level)) {
		DEBUG(std::cout << "LEAF: " << pi << " with t = " << t << std::endl;)
		return processLeaf(t, level, backtrackLevel, completed, groupK, groupL);
	}
	
	const Partition& sigma = *((*pIt).first);
	const RefinementPtr& ref = (*pIt).second;
	++pIt;
	
	uint s = ref->alternatives();
	const bool isBacktrack = ref->type() == Backtrack;
	const bool isGroup = ref->type() == Group;
	const PERM* tForRefinement = &t;
	
	if (isGroup) {
		GroupRefinement<PERM,TRANS>* gref = static_cast<GroupRefinement<PERM,TRANS>*>(ref.get());
		if (this->m_bsgs2 && gref->bsgs() == *this->m_bsgs2) {
			tForRefinement = t2;
		}
	}
	
	ref->sort(*this->m_sorter, &pi);
	typedef typename Refinement<PERM>::RefinementPtrIterator RefIt;
	for (RefIt rIt = ref->backtrackBegin(); rIt != ref->backtrackEnd(); ++rIt) {
		if (isBacktrack && s < groupK.U[backtrackLevel].size()) {
			DEBUG(std::cout << "PRUNE the rest:  s=" << s << " < " << groupK.U[backtrackLevel].size() << std::endl;)
			this->m_statNodesPrunedCosetMinimality += s;
			break;
		}
		
		--s;
		RefinementPtr ref2 = *rIt;
		
		const uint oldFixPointsSize = pi.fixPointsSize();
		DEBUG(std::cout << "  refinement from " << pi << std::endl;)
		const uint strictRefinement = ref2->apply2(pi, *tForRefinement);
		DEBUG(std::cout << "  to " << pi << " with " << strictRefinement << std::endl;)
		DEBUG(for(uint jj=0; jj<level; ++jj) std::cout << " ";)
		DEBUG(std::cout << "NODE " << sigma << "  ~~~>  " << pi << std::endl;)
		/*
		for (uint q = 0; q < level; ++q) std::cout << " ";
		std::cout << " " << level << ": " << sigma << " <-> " << pi2 << " from " << pi << std::endl;
		for (uint q = 0; q < level; ++q) std::cout << " ";
		std::cout << " t = " << t << std::endl;
		*/
		if (!strictRefinement) {
			DEBUG(std::cout << "no strict refinement " << sigma << " -- " << pi << std::endl;)
			++this->m_statNodesPrunedChildRestriction;
			continue;
		}
		if (pi.cells() != sigma.cells()) {
			DEBUG(std::cout << "cell number mismatch " << sigma << " -- " << pi << std::endl;)
			ref2->undo(pi, strictRefinement);
			++this->m_statNodesPrunedChildRestriction;
			continue;
		}
		if (pi.fixPointsSize() != sigma.fixPointsSize()) {
			DEBUG(std::cout << "fix point number mismatch " << sigma << " -- " << pi << std::endl;)
			ref2->undo(pi, strictRefinement);
			++this->m_statNodesPrunedChildRestriction;
			continue;
		}
		PERM tG(t);
		PERM* tH = 0;
		if (pi.fixPointsSize() != oldFixPointsSize) {
			if (!updateMappingPermutation(this->m_bsgs, sigma, pi, tG)) {
				DEBUG(std::cout << "no t found " << sigma << " -- " << pi << "; tG = " << tG << std::endl;)
				ref2->undo(pi, strictRefinement);
				++this->m_statNodesPrunedChildRestriction;
				continue;
			}
			if (this->m_bsgs2) {
				tH = new PERM(*t2);
				if (!updateMappingPermutation(*this->m_bsgs2, sigma, pi, *tH)) {
					DEBUG(std::cout << "no t found " << sigma << " -- " << pi << "; tH = " << tH << std::endl;)
					ref2->undo(pi, strictRefinement);
					++this->m_statNodesPrunedChildRestriction;
					continue;
				}
			}
		}
		if (this->m_pruningLevelDCM && isBacktrack) {
			if (pruneDCM(tG, backtrackLevel, groupK, groupL)) {
				++this->m_statNodesPrunedCosetMinimality2;
				ref2->undo(pi, strictRefinement);
				continue;
			}
		}
		uint ret = search(pIt, pi, tG, tH ? tH : t2, level+1, isBacktrack ? (backtrackLevel + 1) : backtrackLevel, completed, groupK, groupL);
		delete tH;
		DEBUG(std::cout << "retract " << strictRefinement << " from " << pi << " to ";)
		ref2->undo(pi, strictRefinement);
		DEBUG(std::cout <<  pi << std::endl;)
		if (BaseSearch<BSGSIN,TRANSRET>::m_stopAfterFirstElement && ret == 0)
			return 0;
		if (ret < level)
			return ret;
	}
	
	completed = std::min(completed, level);
	return level;
}

template<class BSGSIN,class TRANSRET>
bool RBase<BSGSIN,TRANSRET>::updateMappingPermutation(const BSGSIN& bsgs, const Partition& sigma, const Partition& pi, PERM& t2) const {
	typedef std::vector<ulong>::const_iterator FixIt;
	std::vector<ulong>::const_iterator bIt;
	ulong i = 0;
	FixIt fixSigmaIt = sigma.fixPointsBegin();
	const FixIt fixSigmaEndIt = sigma.fixPointsEnd();
	FixIt fixPiIt = pi.fixPointsBegin();
	DEBUG(print_iterable(bsgs.B.begin(), bsgs.B.end(), 1, "B   ");)
	DEBUG(print_iterable(fixSigmaIt, fixSigmaEndIt, 1, "Sigma");)
	for (bIt = bsgs.B.begin(); bIt != bsgs.B.end(); ++bIt, ++i) {
		DEBUG(std::cout << "  base: " << (*bIt)+1 << std::endl;)
		while (fixSigmaIt != fixSigmaEndIt && *fixSigmaIt != *bIt) {
			DEBUG(std::cout << "  skipping " << (*fixSigmaIt)+1 << " for " << (*bIt)+1 << std::endl;)
			++fixSigmaIt;
			++fixPiIt;
		}
		if (fixSigmaIt == fixSigmaEndIt) {
			DEBUG(std::cout << "  no more fix point found for " << (*bIt)+1 << std::endl;)
			return true;
		}
		const ulong alpha = *fixSigmaIt;
		const ulong beta = *fixPiIt;
		if (t2 / alpha != beta) {
			boost::scoped_ptr<PERM> u_beta(bsgs.U[i].at(t2 % beta));
			if (u_beta) {
				//std::cout << "  multiply with " << *u_beta << " for " << alpha+1 << "," << beta+1 << " // base " << bsgs.B[i] + 1<< std::endl;
				t2 ^= *u_beta;
			} else {
				//std::cout << "could not find a u_b with " << (t2 % beta) << " at " << i << "--" << bsgs.B[i] << " -- " << &bsgs << std::endl;
				return false;
			}
		}
		
		++fixSigmaIt;
		++fixPiIt;
	}
	return true;
}

template<class BSGSIN,class TRANSRET>
const std::vector<ulong>& RBase<BSGSIN,TRANSRET>::subgroupBase() const {
	return m_subgroupBase;
}

}
}

#endif // -- RBASE_H_
