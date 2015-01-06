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


#ifndef PERMUTATION_H_
#define PERMUTATION_H_

#include <permlib/common.h>

// for I/O
#include <string>
#include <iostream>
#include <boost/tokenizer.hpp>
#include <sstream>
#include <set>

#include <boost/shared_ptr.hpp>

namespace permlib {

/// Permutation class storing all values explicitly
class Permutation {
public:
	/// typedef for permutation image
	typedef std::vector<unsigned long> perm;
	
	/// constructs identity permutation acting on n elements
    explicit Permutation(unsigned int n);
	/// constructs permutation acting on n elements, given by string in cycle form
    Permutation(unsigned int n, const std::string &cycles);
	/// sort of copy constructor
    explicit Permutation(const perm &p);
	/// copy constructor
    Permutation(const Permutation &p) : m_perm(p.m_perm), m_isIdentity(p.m_isIdentity) {};

	/// permutation multiplication from the right
    Permutation operator*(const Permutation &p) const;
	/// permutation inplace multiplication from the right
	/**
	 * i.e. THIS := THIS * p
	 */
    Permutation& operator*=(const Permutation &p);
    /// permutation inplace multiplication from the left
	/**
	 * i.e. THIS := p * THIS
	 */
    Permutation& operator^=(const Permutation &p);
	/// permutation inversion
    Permutation operator~() const;
	/// permutation inplace inversion
    Permutation& invertInplace();
	/// equals operator
    bool operator==(const Permutation &p2) const { return m_perm == p2.m_perm; };

	/// lets permutation act on val
    inline ulong operator/(ulong val) const { return at(val); }
	/// lets permutation act on val
    inline ulong at(ulong val) const { return m_perm[val]; }

	/// lets inverse permutation act on val, i.e. compute j such that (this->at(j) == val)
    ulong operator%(ulong val) const;

	/// output in cycle form
    friend std::ostream &operator<< (std::ostream &out, const Permutation &p);

	/// returns true if this permutation is identity
	/**
	 * This is done by checking the image of every point.
	 */
    bool isIdentity() const;
	/// dummy stub for interface compatability with PermutationWord
    inline void flush() {};
	/// number of points this permutation acts on
	inline uint size() const { return m_perm.size(); }
	
	///updates this permutation such that pos is mapped onto val and val onto pos
	void setTransposition(uint pos, uint val);
protected:
	/// defintion of permutation behavior
	perm m_perm;

	/// if set to true, permutation is identity; if set to false then it is not known whether this is identity;
	bool m_isIdentity;

	/// INTERNAL ONLY: constructs an "empty" permutation, i.e. without element mapping
	Permutation(unsigned int n, bool) : m_perm(n), m_isIdentity(false) {}
};


//
//     ----       IMPLEMENTATION
//

inline Permutation::Permutation(unsigned int n) 
	: m_perm(n), m_isIdentity(true) 
{
    for (unsigned int i=0; i<n; ++i)
        m_perm[i] = i;
}

inline Permutation::Permutation(unsigned int n, const std::string & cycles) 
	: m_perm(n), m_isIdentity(false) 
{
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sepCycles(",");
	tokenizer tokens(cycles, sepCycles);

	for (unsigned int i=0; i<n; ++i)
		m_perm[i] = i;

	for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
		std::stringstream ss(*tok_iter);

		ulong first, last, temp;
		ss >> first;
		last = first;

		while (!ss.eof()) {
			ss >> temp;
			m_perm[last-1] = temp-1;
			last = temp;
		}
		m_perm[last-1] = first-1;
	}
}


inline Permutation::Permutation(const perm& p) 
	: m_perm(p), m_isIdentity(false) 
{ }

inline Permutation Permutation::operator*(const Permutation &p) const {
	BOOST_ASSERT(p.m_perm.size() == m_perm.size());

	Permutation res(m_perm.size(), true);
    for (unsigned int i=0; i<m_perm.size(); ++i) {
        res.m_perm[i] = p.m_perm[m_perm[i]];
    }
    return res;
}

inline Permutation& Permutation::operator*=(const Permutation &p) {
	BOOST_ASSERT(p.m_perm.size() == m_perm.size());
	m_isIdentity = false;
	
    for (unsigned int i=0; i<m_perm.size(); ++i) {
        m_perm[i] = p.m_perm[m_perm[i]];
    }
    return *this;
}

inline Permutation& Permutation::operator^=(const Permutation &p) {
	BOOST_ASSERT(p.m_perm.size() == m_perm.size());
	m_isIdentity = false;
	perm tmp(m_perm);

    for (unsigned int i=0; i<m_perm.size(); ++i) {
        m_perm[i] = tmp[p.m_perm[i]];
    }
    return *this;
}

inline Permutation Permutation::operator~() const {
    Permutation res(m_perm.size(), true);
    for (unsigned int i=0; i<m_perm.size(); ++i) {
        res.m_perm[m_perm[i]] = i;
    }
    return res;
}

inline Permutation& Permutation::invertInplace() {
	perm tmp(m_perm);
	for (unsigned int i=0; i<m_perm.size(); ++i) {
		m_perm[tmp[i]] = i;
	}
	return *this;
}

inline ulong Permutation::operator%(ulong val) const {
	for (uint i = 0; i < m_perm.size(); ++i) {
		if (m_perm[i] == val)
			return i;
	}
	// must not happen, we have a permutation!
	BOOST_ASSERT(false);
	return -1;
}

inline bool Permutation::isIdentity() const {
	if (m_isIdentity)
		return true;
	for (unsigned int i=0; i<m_perm.size(); ++i)
		if (at(i) != i)
			return false;
	return true;
}

inline void Permutation::setTransposition(uint pos, uint val) {
	BOOST_ASSERT(pos < m_perm.size());
	BOOST_ASSERT(val < m_perm.size());
	
	m_perm[pos] = val;
	m_perm[val] = pos;
}

inline std::ostream& operator<<(std::ostream& out, const Permutation& p) {
    std::set<ulong> worked;
    bool output = false;
    for (ulong x=0; x<p.m_perm.size(); ++x) {
        unsigned long px, startX;
        startX = x;
        px = p.m_perm[x];
        if (worked.count(x) || x == px) {
            continue;
        }

        worked.insert(x);
        out << "(" << (x+1) << ",";
        while (p.m_perm[px] != startX) {
            out << (px+1) << ",";
            worked.insert(px);
            px = p.m_perm[px];
        }
        worked.insert(px);
        out << (px+1);
        out << ")";
        output = true;
    }
    if (!output)
    	out << "()";
    return out;
}

}

#endif // -- PERMUTATION_H_
