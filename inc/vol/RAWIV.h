/*
  Copyright 2011 The University of Texas at Austin

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

// FFTW-volume RAWIV I/O. As of the libcvc-backed rewrite, the implementation
// in RAWIV.cpp delegates RAWIV header construction, endian handling, and
// per-voxel packing to cvc::volume / cvc::volume_file_info. The byte-swap
// macros and RAWIVHeader struct that previously lived in this header are no
// longer needed and have been removed; the public function signatures below
// are unchanged so existing F2Dock callers compile without modification.

#ifndef __RAWIV_H__
#define __RAWIV_H__

#include "../fft-utils/fftw3.h"
#include "../fft-utils/fftwPrecision.h"

int writeGrid(FFTW_complex *scGrid, FFTW_DATA_TYPE *elecGrid, int n,
              double xCenter, double yCenter, double zCenter, double scale,
              char *fileName, char *fileNameSCRe, char *fileNameSCIm,
              char *fileNameElecRe);
int readRAWIVHeader(int *xDim, int *yDim, int *zDim, double *xCenter,
                    double *yCenter, double *zCenter, double *scale,
                    char *fileName);
int readShapeCompGrid(FFTW_complex **scGrid, int *n, double *xCenter,
                      double *yCenter, double *zCenter, double *scale,
                      char *fileNameRe, char *fileNameIm);
int readShapeCompGrid(FFTW_complex *scGrid, double *xCenter, double *yCenter,
                      double *zCenter, char *fileNameRe, char *fileNameIm);
int readElecGrid(FFTW_DATA_TYPE **elecGrid, int *n, double *xCenter,
                 double *yCenter, double *zCenter, double *scale,
                 char *fileNameRe);
int readElecGrid(FFTW_DATA_TYPE *elecGrid, double *xCenter, double *yCenter,
                 double *zCenter, char *fileNameRe);

#endif
