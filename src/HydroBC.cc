/*
 * HydroBC.cc
 *
 *  Created on: Jan 13, 2012
 *      Author: cferenba
 *
 * Copyright (c) 2012, Los Alamos National Security, LLC.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style open-source
 * license; see top-level LICENSE file for full license text.
 */

#include "HydroBC.hh"

#include "Memory.hh"
#include "Mesh.hh"

using namespace std;


HydroBC::HydroBC(
        Mesh* msh,
        const double2 v,
        const vector<long long>& mbp)
    : mesh(msh), numb(mbp.size()), vfix(v) {

    mapbp = Memory::alloc<long long>(numb);
    copy(mbp.begin(), mbp.end(), mapbp);

    mesh->getPlaneChunks(numb, mapbp, pchbfirst, pchblast);

}


HydroBC::~HydroBC() {}


void HydroBC::applyFixedBC(
        double2* pu,
        double2* pf,
        const long long bfirst,
        const long long blast) {

    #pragma ivdep
    for (long long b = bfirst; b < blast; ++b) {
        long long p = mapbp[b];

        pu[p] = project(pu[p], vfix);
        pf[p] = project(pf[p], vfix);
    }

}

