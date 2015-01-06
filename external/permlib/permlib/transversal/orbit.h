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


#ifndef ORBIT_H_
#define ORBIT_H_

#include <list>

#include <boost/foreach.hpp>

namespace permlib {

/// abstract base class for orbit computation
template<class PERM,class DOMAIN>
class Orbit {
public:
	virtual ~Orbit() {}
	
	/// true iff there exists a transversal element mapping \f$\alpha\f$ to val
	virtual bool contains(const DOMAIN& val) const = 0;

    /// returns one element of the orbit
    virtual const DOMAIN& element() const = 0;
protected:
	/// computes orbit of beta under generators
	/**
	 * @param beta
	 * @param generators
	 * @param a ()-callable structure that defines how a PERM acts on a DOMAIN-element
	 * @param orbitList a list of all orbit elements to be filled by the algorithm
	 */
	template<class Action>
	void orbit(const DOMAIN& beta, const PERMlist &generators, Action a, std::list<DOMAIN>& orbitList);
	
	/// updates an existing orbit of beta after one element has been added
	/**
	 * if this instance of Orbit represents the orbit \f$\beta^{S_0}\f$,
	 * then after call of orbitUpdate it will represent the orbit \f$\beta^{S}\f$ where \f$S = S_0 \cup \{g\}\f$
	 * @param beta
	 * @param generators updated generators, which must include g
	 * @param g new generator which has not been there before
	 * @param a ()-callable structure that defines how a PERM acts on a DOMAIN-element
	 * @param orbitList a list of all orbit elements to be filled by the algorithm
	 */
	template<class Action>
	void orbitUpdate(const DOMAIN& beta, const PERMlist &generators, const PERMptr &g, Action a, std::list<DOMAIN>& orbitList);
	
	/// callback when the orbit algorithm constructs an element alpha_p from alpha and p
	/**
	 * @return true iff alpha_p is a new element that has not been seen before
	 */
	virtual bool foundOrbitElement(const DOMAIN& alpha, const DOMAIN& alpha_p, const PERMptr& p) = 0;
};

template <class PERM,class DOMAIN>
template<class Action>
inline void Orbit<PERM,DOMAIN>::orbit(const DOMAIN& beta, const PERMlist &generators, Action a, std::list<DOMAIN>& orbitList) {
	if (orbitList.empty()) {
		orbitList.push_back(beta);
		foundOrbitElement(beta, beta, PERMptr());
	}
	BOOST_ASSERT( orbitList.size() >= 1 );
	
	DEBUG(std::cout << "orbit of " << beta << std::endl;)
	typename std::list<DOMAIN>::const_iterator it;
	for (it = orbitList.begin(); it != orbitList.end(); ++it) {
		const DOMAIN &alpha = *it;
		BOOST_FOREACH(const PERMptr& p, generators) {
			DOMAIN alpha_p = a(*p, alpha);
			if (foundOrbitElement(alpha, alpha_p, p))
				orbitList.push_back(alpha_p);
		}
	}
}

template <class PERM,class DOMAIN>
template<class Action>
inline void Orbit<PERM,DOMAIN>::orbitUpdate(const DOMAIN& beta, const PERMlist &generators, const PERMptr &g, Action a, std::list<DOMAIN>& orbitList) {
    if (orbitList.empty()) {
		orbitList.push_back(beta);
		foundOrbitElement(beta, beta, PERMptr());
	}
	BOOST_ASSERT( orbitList.size() >= 1 );
	
	DEBUG(std::cout << "orbiUpdate of " << beta << " and " << *g << std::endl;)
	const uint oldSize = orbitList.size();
	// first, compute only ORBIT^g
	typename std::list<DOMAIN>::const_iterator it;
	for (it = orbitList.begin(); it != orbitList.end(); ++it) {
		const DOMAIN &alpha = *it;
		DOMAIN alpha_g = a(*g, alpha);
		if (foundOrbitElement(alpha, alpha_g, g))
			orbitList.push_back(alpha_g);
	}
	
	if (oldSize == orbitList.size())
		return;
	
	orbit(beta, generators, a, orbitList);
}

}

#endif // -- ORBIT_H_
