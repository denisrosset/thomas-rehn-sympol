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


#ifndef STABILIZESPOINTPREDICATE_H_
#define STABILIZESPOINTPREDICATE_H_

#include <boost/foreach.hpp>

namespace permlib {

/// predicate matching points that are stabilized by given permutations
template <class PERM>
class StabilizesPointPredicate : public std::unary_function<ulong, bool> {
public:
	/// constructor
	/**
	 * @param begin begin iterator of permutations (PERM) that the point is checked against
	 * @param end   end   iterator of permutations (PERM) that the point is checked against
	 */
    template<class InputIterator>
    StabilizesPointPredicate(InputIterator begin, InputIterator end) 
		: m_toStabilize(begin, end)
	{ }

	/// evaluate predicate
    bool operator()(const ulong beta) const {
        BOOST_FOREACH(const PERM &p, m_toStabilize) {
            if (p / beta != beta)
                return false;
        }
        return true;
    };
private:
    std::vector<PERM> m_toStabilize;
};

}

#endif // -- STABILIZESPOINTPREDICATE_H_
