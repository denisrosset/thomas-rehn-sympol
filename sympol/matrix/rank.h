// ---------------------------------------------------------------------------
//
// This file is part of SymPol
//
// Copyright (C) 2006-2010  Thomas Rehn <thomas@carmen76.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
// ---------------------------------------------------------------------------

#ifndef SYMPOL_MATRIX_RANK_H_
#define SYMPOL_MATRIX_RANK_H_

#include "algorithm.h"

#include <algorithm>		// for std::min

namespace sympol {
namespace matrix {

/// computes the rank of a matrix by LUP-decomposition/Gaussian eliminiation
template<class Matrix>
class Rank : public Algorithm<Matrix> {
public:
	Rank(Matrix* matrix) : Algorithm<Matrix>(matrix) {}

	using Algorithm<Matrix>::at;
	using Algorithm<Matrix>::m_matrix;

	ulong rank();
};

template<class Matrix>
inline ulong Rank<Matrix>::rank() {
	if (m_matrix->cols() > m_matrix->rows())
		m_matrix->transpose();

	const ulong maxRank = std::min(m_matrix->cols(), m_matrix->rows());
	ulong rank = 0;

	const ulong m = m_matrix->rows();
	const ulong n = m_matrix->cols();
	std::vector<ulong> pi(m);
	for (uint i = 0; i < m; ++i)
		pi[i] = i;

	for (uint k = 0; k < n; ++k) {
		typename Matrix::Type p;
		uint k_prime = 0;
		for (uint i = k; i < m; ++i) {
			if (cmp(abs(at(i,k)), p) > 0) {
				p = abs(at(i,k));
				k_prime = i;
			}
		}
		if (sgn(p) == 0) {
			continue;
		}
		++rank;
		if (rank == maxRank)
			return rank;

		std::swap(pi[k], pi[k_prime]);
		for (uint i = 0; i < n; ++i) {
			std::swap(at(k,i), at(k_prime, i));
		}
		for (uint i = k+1; i < m; ++i) {
			at(i,k) /= at(k,k);
			for (uint j = k+1; j < n; ++j) {
				at(i,j) -= at(i,k) * at(k,j);
			}
		}
	}
	return rank;
}

} // ::matrix
} // ::sympol

#endif // SYMPOL_MATRIX_RANK_H_
