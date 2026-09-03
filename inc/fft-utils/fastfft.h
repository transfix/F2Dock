/*
  Copyright 2011 The University of Texas at Austin

        Authors: Rezaul Alam Chowdhury <shaikat@cs.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of F2Dock.

  F2Dock is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  F2Dock is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
*/

#ifndef CCV_NDFT_H
#define CCV_NDFT_H

#include "fftw3.h"
#include "fftwPrecision.h"
#include "math/SmoothingFunction.h"
#include "sparsefft3-api.h"

void gridding(int M, double *x, double *y, double *z, float *r, char *type,
              FFTW_complex *f, double blobbiness, int n, int m, bool smoothSkin,
              SmoothingFunction *smoothingFunction, FFTW_complex *gHat,
              bool spreadSkin);

void griddingElec(int M, double *x, double *y, double *z, float *r, char *type,
                  FFTW_complex *f, double blobbiness, int n,
                  double elecRadiusInGrids, FFTW_DATA_TYPE *gHat,
                  bool forInPlaceFFT, bool movingMol);

void griddingHbond(int M, double *x, double *y, double *z, float *r,
                   double rExp, FFTW_complex *f, double blobbiness, int n,
                   FFTW_complex *gHat, bool movingMol);

void griddingHydrophobicity(int M, double *x, double *y, double *z, float *r,
                            FFTW_complex *f, double blobbiness, int n,
                            FFTW_complex *gHat, double hydroRadExt,
                            bool pairWise);

void griddingSimpleComplementarity(int M, double *x, double *y, double *z,
                                   float *r, FFTW_complex *f, double blobbiness,
                                   int n, FFTW_complex *gHat,
                                   double simpleRadExt);

bool getCenterFrequencies(int M, double *x, double *y, double *z, float *r,
                          char *type, FFTW_complex *f, double blobbiness,
                          double alpha, int N, int m, FFTW_complex *hHat,
                          bool smoothSkin, SmoothingFunction *smoothingFunction,
                          FFTW_complex *gHat, FFTW_complex *cHat,
                          FFTW_plan gHatPlan, sparse3DFFT_plan gHatSparsePlan,
                          bool griddingDone, bool spreadSkin);

bool getCenterElecFrequencies(int M, double *x, double *y, double *z, float *r,
                              char *type, FFTW_complex *f, double blobbiness,
                              double alpha, int N, double elecRadiusInGrids,
                              FFTW_complex *hHat, FFTW_DATA_TYPE *gHat,
                              FFTW_plan gHatPlan, bool griddingDone,
                              bool movingMol);

bool getCenterHbondFrequencies(int M, double *x, double *y, double *z, float *r,
                               double rExp, FFTW_complex *f, double blobbiness,
                               double alpha, int N, FFTW_complex *hHat,
                               FFTW_complex *gHat, FFTW_complex *cHat,
                               FFTW_plan gHatPlan,
                               sparse3DFFT_plan gHatSparsePlan,
                               bool griddingDone, bool movingMol);

bool getCenterHydrophobicityFrequencies(
    int M, double *x, double *y, double *z, float *r, FFTW_complex *f,
    double blobbiness, double alpha, int N, double hydroRadExt, bool pairWise,
    FFTW_complex *hHat, FFTW_complex *gHat, FFTW_complex *cHat,
    FFTW_plan gHatPlan, sparse3DFFT_plan gHatSparsePlan, bool griddingDone);

bool getCenterSimpleComplementarityFrequencies(
    int M, double *x, double *y, double *z, float *r, FFTW_complex *f,
    double blobbiness, double alpha, int N, double simpleRadExt,
    FFTW_complex *hHat, FFTW_complex *gHat, FFTW_complex *cHat,
    FFTW_plan gHatPlan, sparse3DFFT_plan gHatSparsePlan, bool griddingDone);

#endif
