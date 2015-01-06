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

#ifndef PERMLIB_API_H
#define PERMLIB_API_H

#include <permlib/permutation.h>
#include <permlib/bsgs.h>
#include <permlib/transversal/schreier_tree_transversal.h>
#include <permlib/transversal/orbit_set.h>
#include <permlib/construct/schreier_sims_construction.h>
#include <permlib/change/conjugating_base_change.h>
#include <permlib/search/classic/set_stabilizer_search.h>
#include <permlib/search/classic/set_image_search.h>
#include <permlib/search/orbit_lex_min_search.h>

#include <boost/shared_ptr.hpp>

namespace permlib {

// ---------------------------------------------------------------------
// useful type definitions
//

// need this definition for the ugly PERMlist hack
typedef Permutation PERM;
typedef SchreierTreeTransversal<PERM> TRANSVERSAL;
typedef BSGS<PERM,TRANSVERSAL> PermutationGroup;
typedef OrbitSet<PERM,unsigned long> OrbitAsSet;


// ---------------------------------------------------------------------
// BSGS construction
//

template<class InputIterator>
boost::shared_ptr<PermutationGroup> construct(unsigned long n, InputIterator begin, InputIterator end) {
	SchreierSimsConstruction<PERM, TRANSVERSAL> schreierSims(n);
	boost::shared_ptr<PermutationGroup> group(new PermutationGroup(schreierSims.construct(begin, end)));
	return group;
}


// ---------------------------------------------------------------------
// setwise stabilizer
//

template<class InputIterator>
boost::shared_ptr<PermutationGroup> setStabilizer(const PermutationGroup& group, InputIterator begin, InputIterator end) {
    PermutationGroup copy(group);
	// change the base so that is prefixed by the set
	ConjugatingBaseChange<PERM,TRANSVERSAL,
		RandomBaseTranspose<PERM,TRANSVERSAL> > baseChange(copy);
	baseChange.change(copy, begin, end);
	
	// prepare search without DCM pruning
	classic::SetStabilizerSearch<BSGS<PERM,TRANSVERSAL>, TRANSVERSAL> backtrackSearch(copy, 0);
	backtrackSearch.construct(begin, end);
	
	// start the search
	boost::shared_ptr<PermutationGroup> stabilizer(new PermutationGroup(copy.n));
	backtrackSearch.search(*stabilizer);
	return stabilizer;
}


// ---------------------------------------------------------------------
// set image
//

template<class InputIterator>
boost::shared_ptr<Permutation> setImage(const PermutationGroup& group, InputIterator begin, InputIterator end, InputIterator begin2, InputIterator end2) {
    PermutationGroup copy(group);
	// change the base so that is prefixed by the set
	ConjugatingBaseChange<PERM,TRANSVERSAL,
		RandomBaseTranspose<PERM,TRANSVERSAL> > baseChange(copy);
	baseChange.change(copy, begin, end);
	
	// prepare search without DCM pruning
	classic::SetImageSearch<BSGS<PERM,TRANSVERSAL>, TRANSVERSAL> backtrackSearch(copy, 0);
	backtrackSearch.construct(begin, end, begin2, end2);
	
	// start the search
	return backtrackSearch.searchCosetRepresentative();
}


// ---------------------------------------------------------------------
// orbits
//

std::list<boost::shared_ptr<OrbitAsSet> > orbits(const PermutationGroup& group) {
	std::list<boost::shared_ptr<OrbitAsSet> > orbitList;
	
	for (ulong alpha = 0; alpha < group.n; ++alpha) {
		bool knownElement = false;
		BOOST_FOREACH(const boost::shared_ptr<OrbitAsSet>& orb, orbitList) {
			if (orb->contains(alpha)) {
				knownElement = true;
				break;
			}
		}
		
		if (knownElement)
			continue;
		
		boost::shared_ptr<OrbitAsSet> orbit(new OrbitAsSet());
		orbit->orbit(alpha, group.S, Transversal<PERM>::TrivialAction());
		orbitList.push_back(orbit);
	}

	return orbitList;
}


// ---------------------------------------------------------------------
// smallest orbit element
//

dset smallestSetImage(const PermutationGroup& group, const dset& set) {
	OrbitLexMinSearch<PermutationGroup>  orbLexMin(group);
	return orbLexMin.lexMin(set);
}


} // namespace permlib


#endif // PERMLIB_API_H

