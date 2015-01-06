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


#ifndef TRANSVERSAL_H_
#define TRANSVERSAL_H_

#include <permlib/sorter/base_sorter.h>
#include <permlib/transversal/orbit.h>

#include <map>
#include <list>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

namespace permlib {

template <class PERM>
class Transversal;

template <class PERM>
std::ostream &operator<< (std::ostream &out, const Transversal<PERM> &t) {
	out << "{";
	BOOST_FOREACH (boost::shared_ptr<PERM> p, t.m_transversal) {
		if (p)
			out << *p << ", ";
		else
			out << "O, ";
	}
	out << "}";
	return out;
}

/// Transversal base class corresponding to a base element \f$\alpha\f$
template <class PERM>
class Transversal : public Orbit<PERM,ulong> {
public:
	/// constructor
	/**
	 * @param n size of the set the group is working on
	 */
	Transversal(uint n);
	/// virtual destructor
    virtual ~Transversal() {}
	
	/// returns a transversal element \f$u\f$ such that \f$\alpha^u\f$ equals val
    virtual PERM* at(ulong val) const = 0;
    
	/// true if Schreier generator constructed from x and the transversal element related to "to" is trivial by defintion
	virtual bool trivialByDefinition(const PERM& x, ulong to) const = 0;
	
	/// true iff there exists a transversal element mapping \f$\alpha\f$ to val
	virtual bool contains(const ulong& val) const;

	/// begin iterator of basic orbit
	std::list<ulong>::const_iterator begin() const { return this->m_orbit.begin(); };
	/// end iterator of basic orbit
    std::list<ulong>::const_iterator end() const { return this->m_orbit.end(); };

	/// size of basic orbit / transversal
    uint size() const { return this->m_orbit.size(); }
	
	/// size of the set the group is working on
    inline uint n() const { return m_n; }

	/// sorts orbit according to order given by list of points
	/**
	 * @param Bbegin begin iterator of point list (ulong) inducing an order
	 * @param Bend   end   iterator of point list (ulong) inducing an order
	 */
    template <class InputIterator>
    void sort(InputIterator Bbegin, InputIterator Bend);
    
	/// true iff orbit is sorted
	inline bool sorted() const { return m_sorted; }
	
	/// action of a PERM on ulong element
	struct TrivialAction {
		/// action
		ulong operator()(const PERM &p, ulong v) const {
			return p / v;
		}
	};
	
	/// computes transversal based on orbit of \f$\alpha\f$ under generators
	/**
	 * @param alpha \f$\alpha\f$
	 * @param generators group generators for the orbit
	 */
	virtual void orbit(ulong alpha, const PERMlist &generators);
	/// updates transversal based on orbit of \f$\alpha\f$ under generators where g is a new generator
	/**
	 * @param alpha \f$\alpha\f$
	 * @param generators group generators for the orbit
	 * @param g new generator that the transversal is updated for
	 */
	virtual void orbitUpdate(ulong alpha, const PERMlist &generators, const PERMptr &g);
	
	/// updates transversal after group generators have been conjugated by g
	/**
	 * @param g permutation to conjugate
	 * @param gInv inverse of g for performance reasons
	 */
	virtual void permute(const PERM& g, const PERM& gInv);
	/// updates transversal after group generators have been exchanged
	/**
	 * @param generatorChange map of old generators to new generators
	 */
	virtual void updateGenerators(const std::map<PERM*,PERMptr>& generatorChange) {}
    
    virtual const ulong& element() const;
	
	/// to stream
    friend std::ostream &operator<< <> (std::ostream &out, const Transversal<PERM> &p);
protected:
	/// size of the set the group is working on
    uint m_n;
	
	/// transversal elements
	std::vector<boost::shared_ptr<PERM> > m_transversal;
	
	/// orbit elements
    std::list<ulong> m_orbit;
	
	/// true if orbit is sorted (according to a previous sort(InputIterator, InputIterator) call
    bool m_sorted;
    
	/// stores that 'p' maps 'from' onto 'to'
    virtual void registerMove(ulong from, ulong to, const PERMptr &p);
	
	virtual bool foundOrbitElement(const ulong& alpha, const ulong& alpha_p, const PERMptr& p);
};

//
//     ----       IMPLEMENTATION
//

template <class PERM>
Transversal<PERM>::Transversal(uint n) 
	: m_n(n), m_transversal(n), m_sorted(false) 
{ }

template <class PERM>
void Transversal<PERM>::orbit(ulong beta, const PERMlist &generators) {
	return Orbit<PERM,ulong>::orbit(beta, generators, TrivialAction(), m_orbit);
}

template <class PERM>
void Transversal<PERM>::orbitUpdate(ulong beta, const PERMlist &generators, const PERMptr &g) {
	return Orbit<PERM,ulong>::orbitUpdate(beta, generators, g, TrivialAction(), m_orbit);
}

template <class PERM>
bool Transversal<PERM>::foundOrbitElement(const ulong& alpha, const ulong& alpha_p, const PERMptr& p) {
	if (!m_transversal[alpha_p]) {
		if (!p) {
			PERMptr identity(new PERM(m_n));
			registerMove(alpha, alpha_p, identity);
		} else {
			registerMove(alpha, alpha_p, p);
		}
		return true;
	}
	return false;
}

template <class PERM>
bool Transversal<PERM>::contains(const ulong& val) const { 
	return m_transversal[val] != 0; 
}

template <class PERM>
void Transversal<PERM>::registerMove(ulong from, ulong to, const PERMptr &p) {
	m_sorted = false; 
}


template <class PERM>
template <class InputIterator>
void Transversal<PERM>::sort(InputIterator Bbegin, InputIterator Bend) {
	this->m_orbit.sort(BaseSorter(m_n, Bbegin, Bend));
	m_sorted = true;
}

template <class PERM>
void Transversal<PERM>::permute(const PERM& g, const PERM& gInv) {
	std::vector<boost::shared_ptr<PERM> > temp(m_n);
	for (ulong i=0; i<m_n; ++i) {
		const ulong j = g / i;
		temp[j] = m_transversal[i];
	}
	std::copy(temp.begin(), temp.end(), m_transversal.begin());
	BOOST_FOREACH(ulong& alpha, this->m_orbit) {
		alpha = g / alpha;
	}
	m_sorted = false;
}
    
template <class PERM>
inline const ulong& Transversal<PERM>::element() const {
    return m_orbit.front();
}

}

#endif // -- TRANSVERSAL_H_
