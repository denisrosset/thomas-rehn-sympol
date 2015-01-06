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


#ifndef SYMMETRIC_GROUP_TRANSVERSAL_H_
#define SYMMETRIC_GROUP_TRANSVERSAL_H_

namespace permlib {

template<class PERM>
struct SymmetricGroup;

/// transversal of a symmetric group
template<class PERM>
class SymmetricGroupTransversal {
	public:
		/// constructs a transversal of a symmetric group
		/**
		 * @param sg group
		 * @param basePos position of the element in the group base that this transversal belongs to
		 */
		SymmetricGroupTransversal(const SymmetricGroup<PERM>* sg, uint basePos) : symmetricGroup(sg), basePos(basePos) {}
		
		/// computes a transversal element on demand if one exists
		PERM* at(ulong val) const {
			for (uint i=0; i<basePos; ++i) {
				if ((symmetricGroup->B)[i] == val)
					return 0;
			}
			
			PERM* p = new PERM(symmetricGroup->B.size());
			p->setTransposition((symmetricGroup->B)[basePos],val);
			return p;
		}

		/// size of basic orbit / transversal
		uint size() const { return symmetricGroup->n - basePos; }
	private:
		const SymmetricGroup<PERM>* symmetricGroup;
		uint basePos;
};

}

#endif // SYMMETRIC_GROUP_TRANSVERSAL_H_

