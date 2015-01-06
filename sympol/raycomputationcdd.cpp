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


#include "raycomputationcdd.h"
#include "polyhedron.h"

extern "C" {
	#include <setoper.h>
	#include <cdd.h>
}

#include <ctime>
#include <cstdio>

using namespace yal;
using namespace sympol;

//
// static members
//
bool RayComputationCDD::ms_bInitialized = false;
const char* RayComputationCDD::ms_chName = "CDD";
LoggerPtr RayComputationCDD::logger(Logger::getLogger("RayCompCDD"));

RayComputationCDD::RayComputationCDD() 
	: m_lrs(new RayComputationLRS())
{
}

bool RayComputationCDD::initialize() {
	if (!RayComputationCDD::ms_bInitialized) {
		m_lrs->initialize();
		dd_set_global_constants();
		
		RayComputationCDD::ms_bInitialized = true;       
		return true;
	}
	
	return true;
}

bool RayComputationCDD::finish() {
	if (!RayComputationCDD::ms_bInitialized) {
		return true;
	}
	m_lrs->finish();
	dd_free_global_constants();
	
	RayComputationCDD::ms_bInitialized = false;
	
	return true;
}


/**
 * computes dual description of data into rays
 */
bool RayComputationCDD::dualDescription(const Polyhedron & data, std::vector<FaceWithDataPtr> & rays) const {
	dd_MatrixPtr matrix;
	if (!fillModelCDD(data, matrix))
		return false;
	
	dd_ErrorType err;
	dd_PolyhedraPtr poly = dd_DDMatrix2Poly(matrix, &err);
	if (err != dd_NoError) {
		dd_FreeMatrix(matrix);
		return false;
	}
	dd_MatrixPtr G = dd_CopyGenerators(poly);
	dd_Amatrix A = G->matrix;
	
	if (data.homogeneous()) {
		// cdd doesn't return the cone apex, so add it manually
		QArrayPtr row(new QArray(data.dimension()));
		mpq_set_ui((*row)[0], 1, 1);
		Face face(data.faceDescription(*row));
		if (face.count() == data.rows()) {
			FaceWithDataPtr fdPtr(new FaceWithData(face,row));
			rays.push_back(fdPtr);
		}
	}

	for (uint i = 0; i < static_cast<uint>(G->rowsize); ++i) {
		QArrayPtr row(new QArray(data.dimension()));
		row->initFromArray(A[i]);
		const Face f = data.faceDescription(*row);
		FaceWithDataPtr fdPtr(new FaceWithData(f, row, data.incidenceNumber(f)));
		rays.push_back(fdPtr);
	}
	
	dd_FreePolyhedra(poly);
	dd_FreeMatrix(matrix);
	dd_FreeMatrix(G);
	
	return true;
}

bool RayComputationCDD::firstVertex(const Polyhedron & data, Face & f, QArray & q, bool requireRay) const {
	dd_PolyhedraPtr poly = dd_CreatePolyhedraData(data.rows(), data.dimension());
	if (!poly)
		return false;
	
	poly->representation = dd_Inequality;
	poly->homogeneous = data.homogeneous() ? dd_TRUE : dd_FALSE;
	
	uint i = 0;
	BOOST_FOREACH(const QArray& row, data.rowPair()) {
		for (uint j = 0; j < data.dimension(); ++j) {
			dd_set(poly->A[i][j], row[j]);
		}
		if (data.isLinearity(row))
			poly->EqualityIndex[i+1] = 1;
		++i;
	}
	
	dd_ConePtr cone = dd_ConeDataLoad(poly);
	dd_DDInit(cone);
	dd_boolean found = dd_FALSE;
	dd_FindInitialRays(cone, &found);
	if (found == dd_TRUE) {
		dd_InitialDataSetup(cone);
		dd_RayPtr ray = cone->FirstRay;
		bool vertexFound = false;
		while (ray) {
			if (!ray->feasible) {
				ray = ray->Next;
				continue;
			}
			
			// cdd reduces the cone dimension if it finds
			// linearities (marked with cone->newcol[i] == 0)
			uint j = 0;
			for (i = 1; i <= static_cast<uint>(cone->d_orig); ++i) {
				if (cone->newcol[i]) {
					mpq_set(q[i-1], ray->Ray[j]);
					++j;
				} else {
					mpq_set_ui(q[i-1], 0, 1);
				}
			}
			
			// the following condition is always fulfilled because q is always a ray
			// by construction of cdd
			if (!requireRay || (requireRay && q.isRay())) {
				YALLOG_DEBUG3(logger, "found first vertex " << q);
				f = data.faceDescription(q);
				vertexFound = true;
				break;
			}
			ray = ray->Next;
		}
		
		dd_FreePolyhedra(poly);
		return vertexFound;
	}

	dd_FreePolyhedra(poly);
	return false;
}

bool RayComputationCDD::determineRedundancies(Polyhedron & data, std::list<FaceWithData>& myRays) const {
	dd_MatrixPtr matrix;
	if (!fillModelCDD(data, matrix))
		return false;
	
	std::list<ulong> redundancies;
	dd_ErrorType err;
	dd_rowset redundantRows = dd_RedundantRows(matrix, &err);
	if (err != dd_NoError) {
		dd_FreeMatrix(matrix);
		return false;
	}
	
	for (uint i = 0; i < static_cast<uint>(set_card(redundantRows)); ++i)
		if (set_member(i+1, redundantRows))
			redundancies.push_back(i);
	
	data.addRedundancies(redundancies);
	
	set_free(redundantRows);
	dd_FreeMatrix(matrix);
	
	return true;
}

double RayComputationCDD::estimate(const sympol::Polyhedron& data, std::list<FaceWithData>& rays) const
{
    return m_lrs->estimate(data, rays);
}


bool RayComputationCDD::fillModelCDD(const Polyhedron & data, dd_MatrixPtr& matrix) const {
	matrix = dd_CreateMatrix(data.rows(), data.dimension());
	if (!matrix)
		return false;
	
	matrix->representation = dd_Inequality;
	matrix->numbtype = dd_GetNumberType("rational");
	
	uint i = 0;
	BOOST_FOREACH(const QArray& row, data.rowPair()) {
		for (uint j = 0; j < data.dimension(); ++j) {
			dd_set(matrix->matrix[i][j], row[j]);
		}
		if (data.isLinearity(row))
			set_addelem(matrix->linset, i+1);
		++i;
	}
	
	return true;
}


