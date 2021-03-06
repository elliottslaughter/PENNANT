/*
 * ExportGold.cc
 *
 *  Created on: Mar 1, 2012
 *      Author: cferenba
 *
 * Copyright (c) 2012, Los Alamos National Security, LLC.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style open-source
 * license; see top-level LICENSE file for full license text.
 */

#include "ExportGold.hh"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <numeric>

#include "Parallel.hh"
#include "Vec2.hh"
#include "Mesh.hh"

using namespace std;


ExportGold::ExportGold(Mesh* m) : mesh(m) {}

ExportGold::~ExportGold() {}


void ExportGold::write(
        const string& basename,
        const int cycle,
        const double time,
        const double* zr,
        const double* ze,
        const double* zp) {

    writeCaseFile(basename);

    sortZones();
    writeGeoFile(basename, cycle, time);

    writeVarFile(basename, "zr", zr);
    writeVarFile(basename, "ze", ze);
    writeVarFile(basename, "zp", zp);

}


void ExportGold::writeCaseFile(
        const string& basename) {

    if (Parallel::mype > 0) return;

    // open file
    const string filename = basename + ".case";
    ofstream ofs(filename.c_str());
    if (!ofs.good()) {
        cerr << "Cannot open file " << filename << " for writing"
             << endl;
        exit(1);
    }

    // write case info
    ofs << "#" << endl;
    ofs << "# Created by PENNANT" << endl;
    ofs << "#" << endl;

    ofs << "FORMAT" << endl;
    ofs << "type: ensight gold" << endl;

    ofs << "GEOMETRY" << endl;
    ofs << "model: " << basename << ".geo" << endl;

    ofs << "VARIABLE" << endl;
    ofs << "scalar per element: zr " << basename << ".zr" << endl;
    ofs << "scalar per element: ze " << basename << ".ze" << endl;
    ofs << "scalar per element: zp " << basename << ".zp" << endl;

    ofs.close();

}


void ExportGold::writeGeoFile(
        const string& basename,
        const int cycle,
        const double time) {
    using Parallel::numpe;
    using Parallel::mype;

    // open file
    ofstream ofs;
    if (mype == 0) {
        const string filename = basename + ".geo";
        ofs.open(filename.c_str());
        if (!ofs.good()) {
            cerr << "Cannot open file " << filename << " for writing"
                 << endl;
            exit(1);
        }
    }

    // write general header
    if (mype == 0) {
        ofs << scientific;
        ofs << "cycle = " << setw(8) << cycle << endl;
        ofs << setprecision(8);
        ofs << "t = " << setw(15) << time << endl;
        ofs << "node id assign" << endl;
        ofs << "element id given" << endl;

        // write header for the one "part" (entire mesh)
        ofs << "part" << endl;
        ofs << setw(10) << 1 << endl;
        ofs << "universe" << endl;
    } // if mype == 0

    // gather node info to PE 0
    const long long nump = mesh->nump;
    const double2* px = mesh->px;

    long long gnump = nump;
    Parallel::globalSum(gnump);
    vector<long long> penump(mype == 0 ? numpe : 0);
    Parallel::gather(nump, &penump[0]);
    vector<long long> peoffset(mype == 0 ? numpe + 1 : 1);
    partial_sum(penump.begin(), penump.end(), &peoffset[1]);
    long long offset;
    Parallel::scatter(&peoffset[0], offset);
    vector<double2> gpx(mype == 0 ? gnump : 0);
    Parallel::gatherv(&px[0], nump, &gpx[0], &penump[0]);

    // write node info
    if (mype == 0) {
        ofs << "coordinates" << endl;
        ofs << setw(10) << gnump << endl;
        ofs << setprecision(5);
        for (long long p = 0; p < gnump; ++p)
            ofs << setw(12) << gpx[p].x << endl;
        for (long long p = 0; p < gnump; ++p)
            ofs << setw(12) << gpx[p].y << endl;
        // Ensight expects z-coordinates, so write 0 for those
        for (long long p = 0; p < gnump; ++p)
            ofs << setw(12) << 0. << endl;
    } // if mype

    const int* znump = mesh->znump;
    const long long* mapsp1 = mesh->mapsp1;

    const long long ntris = tris.size();
    const long long nquads = quads.size();
    const long long nothers = others.size();

    if (mype == 0) {
        pentris.resize(numpe);
        penquads.resize(numpe);
        penothers.resize(numpe);
    }
    Parallel::gather(ntris, &pentris[0]);
    Parallel::gather(nquads, &penquads[0]);
    Parallel::gather(nothers, &penothers[0]);

    gntris = accumulate(pentris.begin(), pentris.end(), 0);
    gnquads = accumulate(penquads.begin(), penquads.end(), 0);
    gnothers = accumulate(penothers.begin(), penothers.end(), 0);

    vector<long long> pesizes(mype == 0 ? numpe : 0);

    // gather triangle info to PE 0
    vector<long long> trip(3 * ntris);
    vector<long long> gtris(gntris), gtrip(3 * gntris);
    Parallel::gatherv(&tris[0], ntris, &gtris[0], &pentris[0]);
    if (mype == 0) {
        for (long long pe = 0; pe < numpe; ++pe)
            pesizes[pe] = pentris[pe] * 3;
    }
    for (long long t = 0; t < ntris; ++t) {
        long long z = tris[t];
        long long sbase = mapzs[z];
        for (long long i = 0; i < 3; ++i) {
            trip[t * 3 + i] = mapsp1[sbase + i] + offset;
        }
    }
    Parallel::gatherv(&trip[0], 3 * ntris, &gtrip[0], &pesizes[0]);

    // write triangles
    if (mype == 0 && gntris > 0) {
        ofs << "tria3" << endl;
        ofs << setw(10) << gntris << endl;
        for (long long t = 0; t < gntris; ++t)
            ofs << setw(10) << gtris[t] + 1 << endl;
        for (long long t = 0; t < gntris; ++t) {
            for (long long i = 0; i < 3; ++i)
                ofs << setw(10) << gtrip[t * 3 + i] + 1;
            ofs << endl;
        }
    } // if mype == 0 ...

    // gather quad info to PE 0
    vector<long long> quadp(4 * nquads);
    vector<long long> gquads(gnquads), gquadp(4 * gnquads);
    Parallel::gatherv(&quads[0], nquads, &gquads[0], &penquads[0]);
    if (mype == 0) {
        for (long long pe = 0; pe < numpe; ++pe)
            pesizes[pe] = penquads[pe] * 4;
    }
    for (long long q = 0; q < nquads; ++q) {
        long long z = quads[q];
        long long sbase = mapzs[z];
        for (long long i = 0; i < 4; ++i) {
            quadp[q * 4 + i] = mapsp1[sbase + i] + offset;
        }
    }
    Parallel::gatherv(&quadp[0], 4 * nquads, &gquadp[0], &pesizes[0]);

    // write quads
    if (mype == 0 && gnquads > 0) {
        ofs << "quad4" << endl;
        ofs << setw(10) << gnquads << endl;
        for (long long q = 0; q < gnquads; ++q)
            ofs << setw(10) << gquads[q] + 1 << endl;
        for (long long q = 0; q < gnquads; ++q) {
            for (long long i = 0; i < 4; ++i)
                ofs << setw(10) << gquadp[q * 4 + i] + 1;
            ofs << endl;
        }
    } // if mype == 0 ...

    // gather other info to PE 0
    vector<long long> othernump(nothers), otherp;
    vector<long long> gothers(gnothers), gothernump(gnothers);
    Parallel::gatherv(&others[0], nothers, &gothers[0], &penothers[0]);
    for (long long n = 0; n < nothers; ++n) {
        long long z = others[n];
        long long sbase = mapzs[z];
        othernump[n] = znump[z];
        for (long long i = 0; i < znump[z]; ++i) {
            otherp.push_back(mapsp1[sbase + i] + offset);
        }
    }
    Parallel::gatherv(&othernump[0], nothers, &gothernump[0], &penothers[0]);
    long long size = otherp.size();
    Parallel::gather(size, &pesizes[0]);
    long long gsize = accumulate(pesizes.begin(), pesizes.end(), 0);
    vector<long long> gotherp(gsize);
    Parallel::gatherv(&otherp[0], size, &gotherp[0], &pesizes[0]);

    // write others
    if (mype == 0 && gnothers > 0) {
        ofs << "nsided" << endl;
        ofs << setw(10) << gnothers << endl;
        for (long long n = 0; n < gnothers; ++n)
            ofs << setw(10) << gothers[n] + 1 << endl;
        for (long long n = 0; n < gnothers; ++n)
            ofs << setw(10) << gothernump[n] << endl;
        long long gp = 0;
        for (long long n = 0; n < gnothers; ++n) {
            for (long long i = 0; i < gothernump[n]; ++i)
                ofs << setw(10) << gotherp[gp + i] + 1;
            ofs << endl;
            gp += gothernump[n];
        }
    } // if mype == 0 ...

    if (mype == 0) ofs.close();

}


void ExportGold::writeVarFile(
        const string& basename,
        const string& varname,
        const double* var) {
    using Parallel::mype;

    // open file
    ofstream ofs;
    if (mype == 0) {
        const string filename = basename + "." + varname;
        ofs.open(filename.c_str());
        if (!ofs.good()) {
            cerr << "Cannot open file " << filename << " for writing"
                 << endl;
            exit(1);
        }
    } // if mype == 0

    // write header
    if (mype == 0) {
        ofs << scientific << setprecision(5);
        ofs << varname << endl;
        ofs << "part" << endl;
        ofs << setw(10) << 1 << endl;
    } // if mype == 0

    long long ntris = tris.size();
    long long nquads = quads.size();
    long long nothers = others.size();

    // gather values on triangles to PE 0
    vector<double> tvar(ntris), gtvar(gntris);
    for (long long t = 0; t < ntris; ++t) {
        tvar[t] = var[tris[t]];
    }
    Parallel::gatherv(&tvar[0], ntris, &gtvar[0], &pentris[0]);

    // write values on triangles
    if (mype == 0 && gntris > 0) {
        ofs << "tria3" << endl;
        for (long long t = 0; t < gntris; ++t) {
            ofs << setw(12) << gtvar[t] << endl;
        }
    } // if mype == 0 ...

    // gather values on quads to PE 0
    vector<double> qvar(nquads), gqvar(gnquads);
    for (long long q = 0; q < nquads; ++q) {
        qvar[q] = var[quads[q]];
    }
    Parallel::gatherv(&qvar[0], nquads, &gqvar[0], &penquads[0]);

    // write values on quads
    if (mype == 0 && gnquads > 0) {
        ofs << "quad4" << endl;
        for (long long q = 0; q < gnquads; ++q) {
            ofs << setw(12) << gqvar[q] << endl;
        }
    } // if mype == 0 ...

    // gather values on others to PE 0
    vector<double> ovar(nothers), govar(gnothers);
    for (long long n = 0; n < nothers; ++n) {
        ovar[n] = var[others[n]];
    }
    Parallel::gatherv(&ovar[0], nothers, &govar[0], &penothers[0]);

    // write values on others
    if (mype == 0 && gnothers > 0) {
        ofs << "nsided" << endl;
        for (long long n = 0; n < gnothers; ++n) {
            ofs << setw(12) << govar[n] << endl;
        }
    } // if mype == 0 ...

    if (mype == 0) ofs.close();

}


void ExportGold::sortZones() {

    const long long numz = mesh->numz;
    const int* znump = mesh->znump;

    mapzs.resize(numz);

    // sort zones by size, create an inverse map
    long long scount = 0;
    for (long long z = 0; z < numz; ++z) {
        long long zsize = znump[z];
        if (zsize == 3)
            tris.push_back(z);
        else if (zsize == 4)
            quads.push_back(z);
        else // zsize > 4
            others.push_back(z);
        mapzs[z] = scount;
        scount += zsize;
    } // for z

}

