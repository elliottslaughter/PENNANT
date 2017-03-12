/*
 * GenMesh.hh
 *
 *  Created on: Jun 4, 2013
 *      Author: cferenba
 *
 * Copyright (c) 2013, Los Alamos National Security, LLC.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style open-source
 * license; see top-level LICENSE file for full license text.
 */

#ifndef GENMESH_HH_
#define GENMESH_HH_

#include <string>
#include <vector>
#include "Vec2.hh"

// forward declarations
class InputFile;


class GenMesh {
public:

    std::string meshtype;       // generated mesh type
    long long gnzx, gnzy;             // global number of zones, in x and y
                                // directions
    double lenx, leny;          // length of mesh sides, in x and y
                                // directions
    long long numpex, numpey;         // number of PEs to use, in x and y
                                // directions
    long long mypex, mypey;           // my PE index, in x and y directions
    long long nzx, nzy;               // (local) number of zones, in x and y
                                // directions
    long long zxoffset, zyoffset;     // offsets of local zone array into
                                // global, in x and y directions

    GenMesh(const InputFile* inp);
    ~GenMesh();

    void generate(
            std::vector<double2>& pointpos,
            std::vector<long long>& zonestart,
            std::vector<long long>& zonesize,
            std::vector<long long>& zonepoints,
            std::vector<long long>& slavemstrpes,
            std::vector<long long>& slavemstrcounts,
            std::vector<long long>& slavepoints,
            std::vector<long long>& masterslvpes,
            std::vector<long long>& masterslvcounts,
            std::vector<long long>& masterpoints);

    void generateRect(
            std::vector<double2>& pointpos,
            std::vector<long long>& zonestart,
            std::vector<long long>& zonesize,
            std::vector<long long>& zonepoints,
            std::vector<long long>& slavemstrpes,
            std::vector<long long>& slavemstrcounts,
            std::vector<long long>& slavepoints,
            std::vector<long long>& masterslvpes,
            std::vector<long long>& masterslvcounts,
            std::vector<long long>& masterpoints);

    void generatePie(
            std::vector<double2>& pointpos,
            std::vector<long long>& zonestart,
            std::vector<long long>& zonesize,
            std::vector<long long>& zonepoints,
            std::vector<long long>& slavemstrpes,
            std::vector<long long>& slavemstrcounts,
            std::vector<long long>& slavepoints,
            std::vector<long long>& masterslvpes,
            std::vector<long long>& masterslvcounts,
            std::vector<long long>& masterpoints);

    void generateHex(
            std::vector<double2>& pointpos,
            std::vector<long long>& zonestart,
            std::vector<long long>& zonesize,
            std::vector<long long>& zonepoints,
            std::vector<long long>& slavemstrpes,
            std::vector<long long>& slavemstrcounts,
            std::vector<long long>& slavepoints,
            std::vector<long long>& masterslvpes,
            std::vector<long long>& masterslvcounts,
            std::vector<long long>& masterpoints);

    void calcNumPE();

}; // class GenMesh


#endif /* GENMESH_HH_ */
