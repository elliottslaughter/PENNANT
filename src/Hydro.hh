/*
 * Hydro.hh
 *
 *  Created on: Dec 22, 2011
 *      Author: cferenba
 *
 * Copyright (c) 2012, Los Alamos National Security, LLC.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style open-source
 * license; see top-level LICENSE file for full license text.
 */

#ifndef HYDRO_HH_
#define HYDRO_HH_

#include <string>
#include <vector>

#include "Vec2.hh"

// forward declarations
class InputFile;
class Mesh;
class PolyGas;
class TTS;
class QCS;
class HydroBC;


class Hydro {
public:

    // associated mesh object
    Mesh* mesh;

    // children of this object
    PolyGas* pgas;
    TTS* tts;
    QCS* qcs;
    std::vector<HydroBC*> bcs;

    double cfl;                 // Courant number, limits timestep
    double cflv;                // volume change limit for timestep
    double rinit;               // initial density for main mesh
    double einit;               // initial energy for main mesh
    double rinitsub;            // initial density in subregion
    double einitsub;            // initial energy in subregion
    double uinitradial;         // initial velocity in radial direction
    std::vector<double> bcx;    // x values of x-plane fixed boundaries
    std::vector<double> bcy;    // y values of y-plane fixed boundaries

    double dtrec;               // maximum timestep for hydro
    char msgdtrec[80];          // message:  reason for dtrec

    double2* pu;       // point velocity
    double2* pu0;      // point velocity, start of cycle
    double2* pap;      // point acceleration
    double2* pf;       // point force
    double* pmaswt;    // point mass, weighted by 1/r
    double* cmaswt;    // corner contribution to pmaswt

    double* zm;        // zone mass
    double* zr;        // zone density
    double* zrp;       // zone density, middle of cycle
    double* ze;        // zone specific internal energy
                       // (energy per unit mass)
    double* zetot;     // zone total internal energy
    double* zw;        // zone work done in cycle
    double* zwrate;    // zone work rate
    double* zp;        // zone pressure
    double* zss;       // zone sound speed
    double* zdu;       // zone velocity difference

    double2* sfp;      // side force from pressure
    double2* sfq;      // side force from artificial visc.
    double2* sft;      // side force from tts
    double2* cftot;    // corner force, total from all sources

    Hydro(const InputFile* inp, Mesh* m);
    ~Hydro();

    void init();

    void initRadialVel(
            const double vel,
            const long long pfirst,
            const long long plast);

    void doCycle(const double dt);

    void advPosHalf(
            const double2* px0,
            const double2* pu0,
            const double dt,
            double2* pxp,
            const long long pfirst,
            const long long plast);

    void advPosFull(
            const double2* px0,
            const double2* pu0,
            const double2* pa,
            const double dt,
            double2* px,
            double2* pu,
            const long long pfirst,
            const long long plast);

    void calcCrnrMass(
            const double* zr,
            const double* zarea,
            const double* smf,
            double* cmaswt,
            const long long sfirst,
            const long long slast);

    void sumCrnrForce(
            const double2* sf,
            const double2* sf2,
            const double2* sf3,
            double2* cftot,
            const long long sfirst,
            const long long slast);

    void calcAccel(
            const double2* pf,
            const double* pmass,
            double2* pa,
            const long long pfirst,
            const long long plast);

    void calcRho(
            const double* zm,
            const double* zvol,
            double* zr,
            const long long zfirst,
            const long long zlast);

    void calcWork(
            const double2* sf,
            const double2* sf2,
            const double2* pu0,
            const double2* pu,
            const double2* px0,
            const double dt,
            double* zw,
            double* zetot,
            const long long sfirst,
            const long long slast);

    void calcWorkRate(
            const double* zvol0,
            const double* zvol,
            const double* zw,
            const double* zp,
            const double dt,
            double* zwrate,
            const long long zfirst,
            const long long zlast);

    void calcEnergy(
            const double* zetot,
            const double* zm,
            double* ze,
            const long long zfirst,
            const long long zlast);

    void sumEnergy(
            const double* zetot,
            const double* zarea,
            const double* zvol,
            const double* zm,
            const double* smf,
            const double2* px,
            const double2* pu,
            double& ei,
            double& ek,
            const long long zfirst,
            const long long zlast,
            const long long sfirst,
            const long long slast);

    void calcDtCourant(
            const double* zdl,
            double& dtrec,
            char* msgdtrec,
            const long long zfirst,
            const long long zlast);

    void calcDtVolume(
            const double* zvol,
            const double* zvol0,
            const double dtlast,
            double& dtrec,
            char* msgdtrec,
            const long long zfirst,
            const long long zlast);

    void calcDtHydro(
            const double* zdl,
            const double* zvol,
            const double* zvol0,
            const double dtlast,
            const long long zfirst,
            const long long zlast);

    void getDtHydro(
            double& dtnew,
            std::string& msgdtnew);

    void resetDtHydro();

    void writeEnergyCheck();

}; // class Hydro



#endif /* HYDRO_HH_ */
