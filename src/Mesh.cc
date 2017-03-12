/*
 * Mesh.cc
 *
 *  Created on: Jan 5, 2012
 *      Author: cferenba
 *
 * Copyright (c) 2012, Los Alamos National Security, LLC.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style open-source
 * license; see top-level LICENSE file for full license text.
 */

#include "Mesh.hh"

#include <stdint.h>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "Vec2.hh"
#include "Memory.hh"
#include "Parallel.hh"
#include "InputFile.hh"
#include "GenMesh.hh"
#include "WriteXY.hh"
#include "ExportGold.hh"

using namespace std;


Mesh::Mesh(const InputFile* inp) :
    gmesh(NULL), egold(NULL), wxy(NULL) {

    using Parallel::mype;

    chunksize = inp->getInt("chunksize", 0);
    if (chunksize < 0) {
        if (mype == 0)
            cerr << "Error: bad chunksize " << chunksize << endl;
        exit(1);
    }

    subregion = inp->getDoubleList("subregion", vector<double>());
    if (subregion.size() != 0 && subregion.size() != 4) {
        if (mype == 0)
            cerr << "Error:  subregion must have 4 entries" << endl;
        exit(1);
    }

    writexy = inp->getInt("writexy", 0);
    writegold = inp->getInt("writegold", 0);

    gmesh = new GenMesh(inp);
    wxy = new WriteXY(this);
    egold = new ExportGold(this);

    init();
}


Mesh::~Mesh() {
    delete gmesh;
    delete wxy;
    delete egold;
}


void Mesh::init() {

    // generate mesh
    vector<double2> nodepos;
    vector<long long> cellstart, cellsize, cellnodes;
    vector<long long> slavemstrpes, slavemstrcounts, slavepoints;
    vector<long long> masterslvpes, masterslvcounts, masterpoints;
    gmesh->generate(nodepos, cellstart, cellsize, cellnodes,
            slavemstrpes, slavemstrcounts, slavepoints,
            masterslvpes, masterslvcounts, masterpoints);

    nump = nodepos.size();
    numz = cellstart.size();
    nums = cellnodes.size();
    numc = nums;

    // copy cell sizes to mesh
    znump = Memory::alloc<int>(numz);
    copy(cellsize.begin(), cellsize.end(), znump);

    // populate maps:
    // use the cell* arrays to populate the side maps
    initSides(cellstart, cellsize, cellnodes);
    // release memory from cell* arrays
    cellstart.resize(0);
    cellsize.resize(0);
    cellnodes.resize(0);
    // now populate edge maps using side maps
    initEdges();

    // populate chunk information
    initChunks();

    // create inverse map for corner-to-point gathers
    initInvMap();

    // calculate parallel data structures
    initParallel(slavemstrpes, slavemstrcounts, slavepoints,
            masterslvpes, masterslvcounts, masterpoints);
    // release memory from parallel-related arrays
    slavemstrpes.resize(0);
    slavemstrcounts.resize(0);
    slavepoints.resize(0);
    masterslvpes.resize(0);
    masterslvcounts.resize(0);
    masterpoints.resize(0);

    // write mesh statistics
    writeStats();

    // allocate remaining arrays
    px = Memory::alloc<double2>(nump);
    ex = Memory::alloc<double2>(nume);
    zx = Memory::alloc<double2>(numz);
    px0 = Memory::alloc<double2>(nump);
    pxp = Memory::alloc<double2>(nump);
    exp = Memory::alloc<double2>(nume);
    zxp = Memory::alloc<double2>(numz);
    sarea = Memory::alloc<double>(nums);
    svol = Memory::alloc<double>(nums);
    zarea = Memory::alloc<double>(numz);
    zvol = Memory::alloc<double>(numz);
    sareap = Memory::alloc<double>(nums);
    svolp = Memory::alloc<double>(nums);
    zareap = Memory::alloc<double>(numz);
    zvolp = Memory::alloc<double>(numz);
    zvol0 = Memory::alloc<double>(numz);
    ssurfp = Memory::alloc<double2>(nums);
    elen = Memory::alloc<double>(nume);
    zdl = Memory::alloc<double>(numz);
    smf = Memory::alloc<double>(nums);

    // do a few initial calculations
    #pragma omp parallel for schedule(static)
    for (long long pch = 0; pch < numpch; ++pch) {
        long long pfirst = pchpfirst[pch];
        long long plast = pchplast[pch];
        // copy nodepos into px, distributed across threads
        for (long long p = pfirst; p < plast; ++p)
            px[p] = nodepos[p];

    }

    numsbad = 0;
    #pragma omp parallel for schedule(static)
    for (long long sch = 0; sch < numsch; ++sch) {
        long long sfirst = schsfirst[sch];
        long long slast = schslast[sch];
        calcCtrs(px, ex, zx, sfirst, slast);
        calcVols(px, zx, sarea, svol, zarea, zvol, sfirst, slast);
        calcSideFracs(sarea, zarea, smf, sfirst, slast);
    }
    checkBadSides();

}


void Mesh::initSides(
        const vector<long long>& cellstart,
        const vector<long long>& cellsize,
        const vector<long long>& cellnodes) {

    mapsp1 = Memory::alloc<long long>(nums);
    mapsp2 = Memory::alloc<long long>(nums);
    mapsz  = Memory::alloc<long long>(nums);
    mapss3 = Memory::alloc<long long>(nums);
    mapss4 = Memory::alloc<long long>(nums);

    for (long long z = 0; z < numz; ++z) {
        long long sbase = cellstart[z];
        long long size = cellsize[z];
        for (long long n = 0; n < size; ++n) {
            long long s = sbase + n;
            long long snext = sbase + (n + 1 == size ? 0 : n + 1);
            long long slast = sbase + (n == 0 ? size : n) - 1;
            mapsz[s] = z;
            mapsp1[s] = cellnodes[s];
            mapsp2[s] = cellnodes[snext];
            mapss3[s] = slast;
            mapss4[s] = snext;
        } // for n
    } // for z

}


void Mesh::initEdges() {

    vector<vector<long long> > edgepp(nump), edgepe(nump);

    mapse = Memory::alloc<long long>(nums);

    long long e = 0;
    for (long long s = 0; s < nums; ++s) {
        long long p1 = min(mapsp1[s], mapsp2[s]);
        long long p2 = max(mapsp1[s], mapsp2[s]);

        vector<long long>& vpp = edgepp[p1];
        vector<long long>& vpe = edgepe[p1];
        long long i = find(vpp.begin(), vpp.end(), p2) - vpp.begin();
        if (i == vpp.size()) {
            // (p, p2) isn't in the edge list - add it
            vpp.push_back(p2);
            vpe.push_back(e);
            ++e;
        }
        mapse[s] = vpe[i];
    }  // for s

    nume = e;

}


void Mesh::initChunks() {

    if (chunksize == 0) chunksize = max(nump, nums);

    // compute side chunks
    // use 'chunksize' for maximum chunksize; decrease as needed
    // to ensure that no zone has its sides split across chunk
    // boundaries
    long long s1, s2 = 0;
    while (s2 < nums) {
        s1 = s2;
        s2 = min(s2 + chunksize, nums);
        while (s2 < nums && mapsz[s2] == mapsz[s2-1])
            --s2;
        schsfirst.push_back(s1);
        schslast.push_back(s2);
        schzfirst.push_back(mapsz[s1]);
        schzlast.push_back(mapsz[s2-1] + 1);
    }
    numsch = schsfirst.size();

    // compute point chunks
    long long p1, p2 = 0;
    while (p2 < nump) {
        p1 = p2;
        p2 = min(p2 + chunksize, nump);
        pchpfirst.push_back(p1);
        pchplast.push_back(p2);
    }
    numpch = pchpfirst.size();

    // compute zone chunks
    long long z1, z2 = 0;
    while (z2 < numz) {
        z1 = z2;
        z2 = min(z2 + chunksize, numz);
        zchzfirst.push_back(z1);
        zchzlast.push_back(z2);
    }
    numzch = zchzfirst.size();

}


void Mesh::initInvMap() {
    mappcfirst = Memory::alloc<long long>(nump);
    mapccnext = Memory::alloc<long long>(nums);

    vector<pair<long long, long long> > pcpair(nums);
    for (long long c = 0; c < numc; ++c)
        pcpair[c] = make_pair(mapsp1[c], c);
    sort(pcpair.begin(), pcpair.end());
    for (long long i = 0; i < numc; ++i) {
        long long p = pcpair[i].first;
        long long pp = pcpair[i+1].first;
        long long pm = pcpair[i-1].first;
        long long c = pcpair[i].second;
        long long cp = pcpair[i+1].second;

        if (i == 0 || p != pm)  mappcfirst[p] = c;
        if (i+1 == numc || p != pp)
            mapccnext[c] = -1;
        else
            mapccnext[c] = cp;
    }

}


void Mesh::initParallel(
        const vector<long long>& slavemstrpes,
        const vector<long long>& slavemstrcounts,
        const vector<long long>& slavepoints,
        const vector<long long>& masterslvpes,
        const vector<long long>& masterslvcounts,
        const vector<long long>& masterpoints) {
    if (Parallel::numpe == 1) return;

    nummstrpe = slavemstrpes.size();
    mapmstrpepe = Memory::alloc<long long>(nummstrpe);
    copy(slavemstrpes.begin(), slavemstrpes.end(), mapmstrpepe);
    mstrpenumslv = Memory::alloc<long long>(nummstrpe);
    copy(slavemstrcounts.begin(), slavemstrcounts.end(), mstrpenumslv);
    mapmstrpeslv1 = Memory::alloc<long long>(nummstrpe);
    long long count = 0;
    for (long long mstrpe = 0; mstrpe < nummstrpe; ++mstrpe) {
        mapmstrpeslv1[mstrpe] = count;
        count += mstrpenumslv[mstrpe];
    }
    numslv = slavepoints.size();
    mapslvp = Memory::alloc<long long>(numslv);
    copy(slavepoints.begin(), slavepoints.end(), mapslvp);

    numslvpe = masterslvpes.size();
    mapslvpepe = Memory::alloc<long long>(numslvpe);
    copy(masterslvpes.begin(), masterslvpes.end(), mapslvpepe);
    slvpenumprx = Memory::alloc<long long>(numslvpe);
    copy(masterslvcounts.begin(), masterslvcounts.end(), slvpenumprx);
    mapslvpeprx1 = Memory::alloc<long long>(numslvpe);
    count = 0;
    for (long long slvpe = 0; slvpe < numslvpe; ++slvpe) {
        mapslvpeprx1[slvpe] = count;
        count += slvpenumprx[slvpe];
    }
    numprx = masterpoints.size();
    mapprxp = Memory::alloc<long long>(numprx);
    copy(masterpoints.begin(), masterpoints.end(), mapprxp);

}


void Mesh::writeStats() {

    int64_t gnump = nump;
    // make sure that boundary points aren't double-counted;
    // only count them if they are masters
    if (Parallel::numpe > 1) gnump -= numslv;
    int64_t gnumz = numz;
    int64_t gnums = nums;
    int64_t gnume = nume;
    long long gnumpch = numpch;
    long long gnumzch = numzch;
    long long gnumsch = numsch;

    Parallel::globalSum(gnump);
    Parallel::globalSum(gnumz);
    Parallel::globalSum(gnums);
    Parallel::globalSum(gnume);
    Parallel::globalSum(gnumpch);
    Parallel::globalSum(gnumzch);
    Parallel::globalSum(gnumsch);

    if (Parallel::mype > 0) return;

    cout << "--- Mesh Information ---" << endl;
    cout << "Points:  " << gnump << endl;
    cout << "Zones:  "  << gnumz << endl;
    cout << "Sides:  "  << gnums << endl;
    cout << "Edges:  "  << gnume << endl;
    cout << "Side chunks:  " << gnumsch << endl;
    cout << "Point chunks:  " << gnumpch << endl;
    cout << "Zone chunks:  " << gnumzch << endl;
    cout << "Chunk size:  " << chunksize << endl;
    cout << "------------------------" << endl;

}


void Mesh::write(
        const string& probname,
        const int cycle,
        const double time,
        const double* zr,
        const double* ze,
        const double* zp) {

    if (writexy) {
        if (Parallel::mype == 0)
            cout << "Writing .xy file..." << endl;
        wxy->write(probname, zr, ze, zp);
    }
    if (writegold) {
        if (Parallel::mype == 0) 
            cout << "Writing gold file..." << endl;
        egold->write(probname, cycle, time, zr, ze, zp);
    }

}


vector<long long> Mesh::getXPlane(const double c) {

    vector<long long> mapbp;
    const double eps = 1.e-12;

    for (long long p = 0; p < nump; ++p) {
        if (fabs(px[p].x - c) < eps) {
            mapbp.push_back(p);
        }
    }
    return mapbp;

}


vector<long long> Mesh::getYPlane(const double c) {

    vector<long long> mapbp;
    const double eps = 1.e-12;

    for (long long p = 0; p < nump; ++p) {
        if (fabs(px[p].y - c) < eps) {
            mapbp.push_back(p);
        }
    }
    return mapbp;

}


void Mesh::getPlaneChunks(
        const long long numb,
        const long long* mapbp,
        vector<long long>& pchbfirst,
        vector<long long>& pchblast) {

    pchbfirst.resize(0);
    pchblast.resize(0);

    // compute boundary point chunks
    // (boundary points contained in each point chunk)
    long long bf, bl = 0;
    for (long long pch = 0; pch < numpch; ++pch) {
         long long pl = pchplast[pch];
         bf = bl;
         bl = lower_bound(&mapbp[bf], &mapbp[numb], pl) - &mapbp[0];
         pchbfirst.push_back(bf);
         pchblast.push_back(bl);
    }

}


void Mesh::calcCtrs(
        const double2* px,
        double2* ex,
        double2* zx,
        const long long sfirst,
        const long long slast) {

    long long zfirst = mapsz[sfirst];
    long long zlast = (slast < nums ? mapsz[slast] : numz);
    fill(&zx[zfirst], &zx[zlast], double2(0., 0.));

    for (long long s = sfirst; s < slast; ++s) {
        long long p1 = mapsp1[s];
        long long p2 = mapsp2[s];
        long long e = mapse[s];
        long long z = mapsz[s];
        ex[e] = 0.5 * (px[p1] + px[p2]);
        zx[z] += px[p1];
    }

    for (long long z = zfirst; z < zlast; ++z) {
        zx[z] /= (double) znump[z];
    }

}


void Mesh::calcVols(
        const double2* px,
        const double2* zx,
        double* sarea,
        double* svol,
        double* zarea,
        double* zvol,
        const long long sfirst,
        const long long slast) {

    long long zfirst = mapsz[sfirst];
    long long zlast = (slast < nums ? mapsz[slast] : numz);
    fill(&zvol[zfirst], &zvol[zlast], 0.);
    fill(&zarea[zfirst], &zarea[zlast], 0.);

    const double third = 1. / 3.;
    long long count = 0;
    for (long long s = sfirst; s < slast; ++s) {
        long long p1 = mapsp1[s];
        long long p2 = mapsp2[s];
        long long z = mapsz[s];

        // compute side volumes, sum to zone
        double sa = 0.5 * cross(px[p2] - px[p1], zx[z] - px[p1]);
        double sv = third * sa * (px[p1].x + px[p2].x + zx[z].x);
        sarea[s] = sa;
        svol[s] = sv;
        zarea[z] += sa;
        zvol[z] += sv;

        // check for negative side volumes
        if (sv <= 0.) count += 1;

    } // for s

    if (count > 0) {
        #pragma omp atomic
        numsbad += count;
    }

}


void Mesh::checkBadSides() {

    // if there were negative side volumes, error exit
    if (numsbad > 0) {
        cerr << "Error: " << numsbad << " negative side volumes" << endl;
        cerr << "Exiting..." << endl;
        exit(1);
    }

}


void Mesh::calcSideFracs(
        const double* sarea,
        const double* zarea,
        double* smf,
        const long long sfirst,
        const long long slast) {

    #pragma ivdep
    for (long long s = sfirst; s < slast; ++s) {
        long long z = mapsz[s];
        smf[s] = sarea[s] / zarea[z];
    }
}


void Mesh::calcSurfVecs(
        const double2* zx,
        const double2* ex,
        double2* ssurf,
        const long long sfirst,
        const long long slast) {

    #pragma ivdep
    for (long long s = sfirst; s < slast; ++s) {
        long long z = mapsz[s];
        long long e = mapse[s];

        ssurf[s] = rotateCCW(ex[e] - zx[z]);

    }

}


void Mesh::calcEdgeLen(
        const double2* px,
        double* elen,
        const long long sfirst,
        const long long slast) {

    for (long long s = sfirst; s < slast; ++s) {
        const long long p1 = mapsp1[s];
        const long long p2 = mapsp2[s];
        const long long e = mapse[s];

        elen[e] = length(px[p2] - px[p1]);

    }
}


void Mesh::calcCharLen(
        const double* sarea,
        double* zdl,
        const long long sfirst,
        const long long slast) {

    long long zfirst = mapsz[sfirst];
    long long zlast = (slast < nums ? mapsz[slast] : numz);
    fill(&zdl[zfirst], &zdl[zlast], 1.e99);

    for (long long s = sfirst; s < slast; ++s) {
        long long z = mapsz[s];
        long long e = mapse[s];

        double area = sarea[s];
        double base = elen[e];
        double fac = (znump[z] == 3 ? 3. : 4.);
        double sdl = fac * area / base;
        zdl[z] = min(zdl[z], sdl);
    }
}


template <typename T>
void Mesh::parallelGather(
        const T* pvar,
        T* prxvar) {
#ifdef USE_MPI
    // This routine gathers slave values for which MYPE owns the masters.
    const long long tagmpi = 100;
    const long long type_size = sizeof(T);
//    std::vector<T> slvvar(numslv);
    T* slvvar = Memory::alloc<T>(numslv);

    // Post receives for incoming messages from slaves.
    // Store results in proxy buffer.
//    vector<MPI_Request> request(numslvpe);
    MPI_Request* request = Memory::alloc<MPI_Request>(numslvpe);
    for (long long slvpe = 0; slvpe < numslvpe; ++slvpe) {
        long long pe = mapslvpepe[slvpe];
        long long nprx = slvpenumprx[slvpe];
        long long prx1 = mapslvpeprx1[slvpe];
        MPI_Irecv(&prxvar[prx1], nprx * type_size, MPI_BYTE,
                pe, tagmpi, MPI_COMM_WORLD, &request[slvpe]);
    }

    // Load slave data buffer from points.
    for (long long slv = 0; slv < numslv; ++slv) {
        long long p = mapslvp[slv];
        slvvar[slv] = pvar[p];
    }

    // Send slave data to master PEs.
    for (long long mstrpe = 0; mstrpe < nummstrpe; ++mstrpe) {
        long long pe = mapmstrpepe[mstrpe];
        long long nslv = mstrpenumslv[mstrpe];
        long long slv1 = mapmstrpeslv1[mstrpe];
        MPI_Send(&slvvar[slv1], nslv * type_size, MPI_BYTE,
                pe, tagmpi, MPI_COMM_WORLD);
    }

    // Wait for all receives to complete.
//    vector<MPI_Status> status(numslvpe);
    MPI_Status* status = Memory::alloc<MPI_Status>(numslvpe);
    long long ierr = MPI_Waitall(numslvpe, &request[0], &status[0]);
    if (ierr != 0) {
        cerr << "Error: parallelGather MPI error " << ierr <<
                " on PE " << Parallel::mype << endl;
        cerr << "Exiting..." << endl;
        exit(1);
    }

    Memory::free(slvvar);
    Memory::free(request);
    Memory::free(status);
#endif
}


template <typename T>
void Mesh::parallelSum(
        T* pvar,
        T* prxvar) {
#ifdef USE_MPI
    // Compute sum of all (proxy/master) sets.
    // Store results in master.
    for (long long prx = 0; prx < numprx; ++prx) {
        long long p = mapprxp[prx];
        pvar[p] += prxvar[prx];
    }

    // Copy updated master data back to proxies.
    for (long long prx = 0; prx < numprx; ++prx) {
        long long p = mapprxp[prx];
        prxvar[prx] = pvar[p];
    }
#endif
}


template <typename T>
void Mesh::parallelScatter(
        T* pvar,
        const T* prxvar) {
#ifdef USE_MPI
    // This routine scatters master values on MYPE to all slave copies
    // owned by other PEs.
    const long long tagmpi = 200;
    const long long type_size = sizeof(T);
//    std::vector<T> slvvar(numslv);
    T* slvvar = Memory::alloc<T>(numslv);

    // Post receives for incoming messages from masters.
    // Store results in slave buffer.
//    vector<MPI_Request> request(nummstrpe);
    MPI_Request* request = Memory::alloc<MPI_Request>(nummstrpe);
    for (long long mstrpe = 0; mstrpe < nummstrpe; ++mstrpe) {
        long long pe = mapmstrpepe[mstrpe];
        long long nslv = mstrpenumslv[mstrpe];
        long long slv1 = mapmstrpeslv1[mstrpe];
        MPI_Irecv(&slvvar[slv1], nslv * type_size, MPI_BYTE,
                pe, tagmpi, MPI_COMM_WORLD,  &request[mstrpe]);
    }

    // Send updated slave data from proxy buffer back to slave PEs.
    for (long long slvpe = 0; slvpe < numslvpe; ++slvpe) {
        long long pe = mapslvpepe[slvpe];
        long long nprx = slvpenumprx[slvpe];
        long long prx1 = mapslvpeprx1[slvpe];
        MPI_Send((void*)&prxvar[prx1], nprx * type_size, MPI_BYTE,
                pe, tagmpi, MPI_COMM_WORLD);
    }

    // Wait for all receives to complete.
//    vector<MPI_Status> status(nummstrpe);
    MPI_Status* status = Memory::alloc<MPI_Status>(nummstrpe);
    long long ierr = MPI_Waitall(nummstrpe, &request[0], &status[0]);
    if (ierr != 0) {
        cerr << "Error: parallelScatter MPI error " << ierr <<
                " on PE " << Parallel::mype << endl;
        cerr << "Exiting..." << endl;
        exit(1);
    }

    // Store slave data from buffer back to points.
    for (long long slv = 0; slv < numslv; ++slv) {
        long long p = mapslvp[slv];
        pvar[p] = slvvar[slv];
    }

    Memory::free(slvvar);
    Memory::free(request);
    Memory::free(status);
#endif
}


template <typename T>
void Mesh::sumAcrossProcs(T* pvar) {
    if (Parallel::numpe == 1) return;
//    std::vector<T> prxvar(numprx);
    T* prxvar = Memory::alloc<T>(numprx);
    parallelGather(pvar, &prxvar[0]);
    parallelSum(pvar, &prxvar[0]);
    parallelScatter(pvar, &prxvar[0]);
    Memory::free(prxvar);
}


template <typename T>
void Mesh::sumOnProc(
        const T* cvar,
        T* pvar) {

    #pragma omp parallel for schedule(static)
    for (long long pch = 0; pch < numpch; ++pch) {
        long long pfirst = pchpfirst[pch];
        long long plast = pchplast[pch];
        for (long long p = pfirst; p < plast; ++p) {
            T x = T();
            for (long long c = mappcfirst[p]; c >= 0; c = mapccnext[c]) {
                x += cvar[c];
            }
            pvar[p] = x;
        }  // for p
    }  // for pch

}


template <>
void Mesh::sumToPoints(
        const double* cvar,
        double* pvar) {

    sumOnProc(cvar, pvar);
    if (Parallel::numpe > 1)
        sumAcrossProcs(pvar);

}


template <>
void Mesh::sumToPoints(
        const double2* cvar,
        double2* pvar) {

    sumOnProc(cvar, pvar);
    if (Parallel::numpe > 1)
        sumAcrossProcs(pvar);

}

