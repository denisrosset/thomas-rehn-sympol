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


#ifndef RANDOMSCHREIERSIMSCONSTRUCTION_H_
#define RANDOMSCHREIERSIMSCONSTRUCTION_H_

#include <permlib/construct/base_construction.h>
#include <permlib/bsgs.h>

#include <boost/foreach.hpp>

namespace permlib {

/// BSGS construction with Random Schreier-Sims algorithm
/**
 * Randomized algorithm for BSGS construction. If order is known it is Las Vegas type, otherwise Monte Carlo 
 * (it may return an incomplete BSGS).
 */
template <class PERM, class TRANS>
class RandomSchreierSimsConstruction : public BaseConstruction<PERM, TRANS> {
public:
	/// constructor
	/**
	 * @param n cardinality of the set the group is acting on
	 * @param rng a RandomGenerator generating uniformly distributed random group elements of the group that the BSGS is constructed of
	 * @param knownOrder order of the group that the BSGS is constructed of. If non-zero upgrades algorithm to Las Vegas type and the output is guaranteed to be a BSGS.
	 */
	RandomSchreierSimsConstruction(uint n, RandomGenerator<PERM> *rng, ulong knownOrder = 0);

	/// constructs a probable BSGS for group given by generators with no base prescribed
	/** @see construct(ForwardIterator generatorsBegin, ForwardIterator generatorsEnd, InputIterator prescribedBaseBegin, InputIterator prescribedBaseEnd, bool& guaranteedBSGS)
	 */
	template <class ForwardIterator>
	BSGS<PERM, TRANS> construct(ForwardIterator generatorsBegin, ForwardIterator generatorsEnd, bool& guaranteedBSGS) const;

	/// constructs a probable BSGS for group given by generators respecting prescribed base elements
	/**
	 * runs (#{size of initial base} * maxIterationFactor) iterations or until constructed BSGS has knownOrder 
	 * or #minimalConsecutiveSiftingElementCount sift through the constructed BSGS consecutively
	 *
	 * @param generatorsBegin begin iterator of group generators of type PERM
	 * @param generatorsEnd  end iterator of group generators of type PERM
	 * @param prescribedBaseBegin begin iterator of prescribed base of type ulong
	 * @param prescribedBaseEnd  end iterator of prescribed base of type ulong
	 * @param guaranteedBSGS iff true, return object is guaranteed to be a BSGS
	 */
	template <class ForwardIterator, class InputIterator>
	BSGS<PERM, TRANS> construct(ForwardIterator generatorsBegin, ForwardIterator generatorsEnd, InputIterator prescribedBaseBegin, InputIterator prescribedBaseEnd, bool& guaranteedBSGS) const;

	/// number of Schreier generators examined during the last construct call
	mutable uint m_statRandomElementsConsidered;

	/// number of elements that have to sift through constructed BSGS consecutively that it is returned as a probable BSGS
	static const uint minimalConsecutiveSiftingElementCount = 20;
	
	/// factor limiting the number of maximal iterations depeding on the initial base size
	static const uint maxIterationFactor = 10000;
private:
	RandomGenerator<PERM> *m_rng;
	ulong m_knownOrder;
};

//
//     ----       IMPLEMENTATION
//

template <class PERM, class TRANS>
RandomSchreierSimsConstruction<PERM,TRANS>::RandomSchreierSimsConstruction(uint n, RandomGenerator<PERM> *rng, ulong knownOrder) 
	: BaseConstruction<PERM, TRANS>(n), m_statRandomElementsConsidered(0), m_rng(rng), m_knownOrder(knownOrder)
{ }

template <class PERM, class TRANS>
template <class ForwardIterator>
inline BSGS<PERM, TRANS> RandomSchreierSimsConstruction<PERM,TRANS>::construct(ForwardIterator generatorsBegin, ForwardIterator generatorsEnd, bool& guaranteedBSGS) const {
	return construct(generatorsBegin, generatorsEnd, BaseConstruction<PERM,TRANS>::empty, BaseConstruction<PERM,TRANS>::empty, guaranteedBSGS);
}

template <class PERM, class TRANS>
template <class ForwardIterator, class InputIterator>
BSGS<PERM, TRANS> RandomSchreierSimsConstruction<PERM, TRANS>
	::construct(ForwardIterator generatorsBegin, ForwardIterator generatorsEnd, InputIterator prescribedBaseBegin, InputIterator prescribedBaseEnd, bool& guaranteedBSGS) const
{
	const uint &n = BaseConstruction<PERM, TRANS>::m_n;
	BSGS<PERM, TRANS> ret(n);
	std::vector<ulong> &B = ret.B;
	std::vector<TRANS> &U = ret.U;
	std::vector<PERMlist> S;
	setup(generatorsBegin, generatorsEnd, prescribedBaseBegin, prescribedBaseEnd, ret, S);
	
	uint consecutiveSiftingElementCount = minimalConsecutiveSiftingElementCount;
	if (m_knownOrder > 0) 
		// remove consecutive sifting limit if we have the group order as Las Vegas-abort criterion
		consecutiveSiftingElementCount = maxIterationFactor;
	const uint maxIterationCount = B.size() * maxIterationFactor;
	for (uint it = 0; it < maxIterationCount; ++it) {
		bool isProbableBSGS = true;
		for (uint i = 0; i < consecutiveSiftingElementCount && ret.order() != m_knownOrder; ++i) {
			PERM g = m_rng->next();
			++m_statRandomElementsConsidered;
			PERM h(n);
			uint j = ret.sift(g, h);
			if (j < B.size() || !h.isIdentity()) {
				// flush it, because we add it as a generator
				h.flush();
				
				if (j == B.size()) {
					ulong gamma = n+1;
					if (ret.chooseBaseElement(h, gamma)) {
						B.push_back(gamma);
					}
					BOOST_ASSERT(j < B.size());
					S.push_back(PERMlist());
					U.push_back(TRANS(n));
				}
				
				boost::shared_ptr<PERM> hPtr(new PERM(h));
				S[j].insert(S[j].end(), hPtr);

				ret.orbitUpdate(j, S[j], hPtr);
				isProbableBSGS = false;
				break;
			}
		}
		if (isProbableBSGS)
			break;
	}
	
	mergeGenerators(S, ret);
	
	// convinience check of group order
	guaranteedBSGS = ret.order() == m_knownOrder;
	
	return ret;
}

}

#endif // -- RANDOMSCHREIERSIMSCONSTRUCTION_H_
