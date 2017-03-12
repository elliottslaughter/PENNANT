/*
 * QCS.hh
 *
 *  Created on: Feb 21, 2012
 *      Author: cferenba
 *
 * Copyright (c) 2012, Los Alamos National Security, LLC.
 * All rights reserved.
 * Use of this source code is governed by a BSD-style open-source
 * license; see top-level LICENSE file for full license text.
 */

#ifndef QCS_HH_
#define QCS_HH_

#include "Vec2.hh"

// forward declarations
class InputFile;
class Hydro;


class QCS {
public:

    // parent hydro object
    Hydro* hydro;

    double qgamma;                 // gamma coefficient for Q model
    double q1, q2;                 // linear and quadratic coefficients
                                   // for Q model

    QCS(const InputFile* inp, Hydro* h);
    ~QCS();

    void calcForce(
            double2* sf,
            const long long sfirst,
            const long long slast);

    void setCornerDiv(
            double* c0area,
            double* c0div,
            double* c0evol,
            double* c0du,
            double* c0cos,
            const long long sfirst,
            const long long slast);

    void setQCnForce(
            const double* c0div,
            const double* c0du,
            const double* c0evol,
            double2* c0qe,
            const long long sfirst,
            const long long slast);

    void setForce(
            const double* c0area,
            const double2* c0qe,
            double* c0cos,
            double2* sfqq,
            const long long sfirst,
            const long long slast);

    void setVelDiff(
            const long long sfirst,
            const long long slast);

};  // class QCS


#endif /* QCS_HH_ */
