/*
 * HydroBC.hh
 *
 *  Created on: Jan 13, 2012
 *      Author: cferenba
 *
 * Copyright (c) 2012, Los Alamos National Security, LLC.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style open-source
 * license; see top-level LICENSE file for full license text.
 */

#ifndef HYDROBC_HH_
#define HYDROBC_HH_

#include <vector>

#include "Vec2.hh"

// forward declarations
class Mesh;


class HydroBC {
public:

    // associated mesh object
    Mesh* mesh;

    long long numb;                      // number of bdy points
    double2 vfix;                  // vector perp. to fixed plane
    long long* mapbp;                    // map: bdy point -> point
    std::vector<long long> pchbfirst;    // start/stop index for bdy pt chunks
    std::vector<long long> pchblast;

    HydroBC(
            Mesh* msh,
            const double2 v,
            const std::vector<long long>& mbp);

    ~HydroBC();

    void applyFixedBC(
            double2* pu,
            double2* pf,
            const long long bfirst,
            const long long blast);

}; // class HydroBC


#endif /* HYDROBC_HH_ */
