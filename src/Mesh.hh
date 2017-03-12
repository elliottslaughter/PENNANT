/*
 * Mesh.hh
 *
 *  Created on: Jan 5, 2012
 *      Author: cferenba
 *
 * Copyright (c) 2012, Los Alamos National Security, LLC.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style open-source
 * license; see top-level LICENSE file for full license text.
 */

#ifndef MESH_HH_
#define MESH_HH_

#include <string>
#include <vector>

#include "Vec2.hh"

// forward declarations
class InputFile;
class GenMesh;
class WriteXY;
class ExportGold;


class Mesh {
public:

    // children
    GenMesh* gmesh;
    WriteXY* wxy;
    ExportGold* egold;

    // parameters
    int chunksize;                 // max size for processing chunks
    std::vector<double> subregion; // bounding box for a subregion
                                   // if nonempty, should have 4 entries:
                                   // xmin, xmax, ymin, ymax
    bool writexy;                  // flag:  write .xy file?
    bool writegold;                // flag:  write Ensight file?

    // mesh variables
    // (See documentation for more details on the mesh
    //  data structures...)
    long long nump, nume, numz, nums, numc;
                       // number of points, edges, zones,
                       // sides, corners, resp.
    long long numsbad;       // number of bad sides (negative volume)
    long long* mapsp1;       // maps: side -> points 1 and 2
    long long* mapsp2;
    long long* mapsz;        // map: side -> zone
    long long* mapse;        // map: side -> edge
    long long* mapss3;       // map: side -> previous side
    long long* mapss4;       // map: side -> next side

    // point-to-corner inverse map is stored as a linked list...
    long long* mappcfirst;   // map:  point -> first corner
    long long* mapccnext;    // map:  corner -> next corner

    // mpi comm variables
    long long nummstrpe;     // number of messages mype sends to master pes
    long long numslvpe;      // number of messages mype receives from slave pes
    long long numprx;        // number of proxies on mype
    long long numslv;        // number of slaves on mype
    long long* mapslvpepe;   // map: slave pe -> (global) pe
    long long* mapslvpeprx1; // map: slave pe -> first proxy in proxy buffer
    long long* mapprxp;      // map: proxy -> corresponding (master) point
    long long* slvpenumprx;  // number of proxies for each slave pe
    long long* mapmstrpepe;  // map: master pe -> (global) pe
    long long* mstrpenumslv; // number of slaves for each master pe
    long long* mapmstrpeslv1;// map: master pe -> first slave in slave buffer
    long long* mapslvp;      // map: slave -> corresponding (slave) point

    int* znump;        // number of points in zone

    double2* px;       // point coordinates
    double2* ex;       // edge center coordinates
    double2* zx;       // zone center coordinates
    double2* pxp;      // point coords, middle of cycle
    double2* exp;      // edge ctr coords, middle of cycle
    double2* zxp;      // zone ctr coords, middle of cycle
    double2* px0;      // point coords, start of cycle

    double* sarea;     // side area
    double* svol;      // side volume
    double* zarea;     // zone area
    double* zvol;      // zone volume
    double* sareap;    // side area, middle of cycle
    double* svolp;     // side volume, middle of cycle
    double* zareap;    // zone area, middle of cycle
    double* zvolp;     // zone volume, middle of cycle
    double* zvol0;     // zone volume, start of cycle

    double2* ssurfp;   // side surface vector
    double* elen;      // edge length
    double* smf;       // side mass fraction
    double* zdl;       // zone characteristic length

    long long numsch;                    // number of side chunks
    std::vector<long long> schsfirst;    // start/stop index for side chunks
    std::vector<long long> schslast;
    std::vector<long long> schzfirst;    // start/stop index for zone chunks
    std::vector<long long> schzlast;
    long long numpch;                    // number of point chunks
    std::vector<long long> pchpfirst;    // start/stop index for point chunks
    std::vector<long long> pchplast;
    long long numzch;                    // number of zone chunks
    std::vector<long long> zchzfirst;    // start/stop index for zone chunks
    std::vector<long long> zchzlast;

    Mesh(const InputFile* inp);
    ~Mesh();

    void init();

    // populate mapping arrays
    void initSides(
            const std::vector<long long>& cellstart,
            const std::vector<long long>& cellsize,
            const std::vector<long long>& cellnodes);
    void initEdges();

    // populate chunk information
    void initChunks();

    // populate inverse map
    void initInvMap();

    void initParallel(
            const std::vector<long long>& slavemstrpes,
            const std::vector<long long>& slavemstrcounts,
            const std::vector<long long>& slavepoints,
            const std::vector<long long>& masterslvpes,
            const std::vector<long long>& masterslvcounts,
            const std::vector<long long>& masterpoints);

    // write mesh statistics
    void writeStats();

    // write mesh
    void write(
            const std::string& probname,
            const int cycle,
            const double time,
            const double* zr,
            const double* ze,
            const double* zp);

    // find plane with constant x, y value
    std::vector<long long> getXPlane(const double c);
    std::vector<long long> getYPlane(const double c);

    // compute chunks for a given plane
    void getPlaneChunks(
            const long long numb,
            const long long* mapbp,
            std::vector<long long>& pchbfirst,
            std::vector<long long>& pchblast);

    // compute edge, zone centers
    void calcCtrs(
            const double2* px,
            double2* ex,
            double2* zx,
            const long long sfirst,
            const long long slast);

    // compute side, corner, zone volumes
    void calcVols(
            const double2* px,
            const double2* zx,
            double* sarea,
            double* svol,
            double* zarea,
            double* zvol,
            const long long sfirst,
            const long long slast);

    // check to see if previous volume computation had any
    // sides with negative volumes
    void checkBadSides();

    // compute side mass fractions
    void calcSideFracs(
            const double* sarea,
            const double* zarea,
            double* smf,
            const long long sfirst,
            const long long slast);

    // compute surface vectors for median mesh
    void calcSurfVecs(
            const double2* zx,
            const double2* ex,
            double2* ssurf,
            const long long sfirst,
            const long long slast);

    // compute edge lengths
    void calcEdgeLen(
            const double2* px,
            double* elen,
            const long long sfirst,
            const long long slast);

    // compute characteristic lengths
    void calcCharLen(
            const double* sarea,
            double* zdl,
            const long long sfirst,
            const long long slast);

    // sum corner variables to points (double or double2)
    template <typename T>
    void sumToPoints(
            const T* cvar,
            T* pvar);

    // helper routines for sumToPoints
    template <typename T>
    void sumOnProc(
            const T* cvar,
            T* pvar);
    template <typename T>
    void sumAcrossProcs(T* pvar);
    template <typename T>
    void parallelGather(
            const T* pvar,
            T* prxvar);
    template <typename T>
    void parallelSum(
            T* pvar,
            T* prxvar);
    template <typename T>
    void parallelScatter(
            T* pvar,
            const T* prxvar);

}; // class Mesh



#endif /* MESH_HH_ */
