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

#ifndef RANKFFTW_H_
#define RANKFFTW_H_

#include <stdio.h>
#include <stdlib.h>

#if !defined(__APPLE__)
#include <malloc.h>
#endif

#include "fftw3.h"
#include "fftwPrecision.h"
#include "../utils/utils.h"

#define NO_ERROR 0
#define INVALID_PARAMETER -1
#define MEMORY_ALLOCATION_FAILED -2
#define FFTW_FORWARD_FAILED -3
#define FFTW_BACKWARD_FAILED -4
#define FFTW_FORWARD_WISDOM_FAILED -5
#define FFTW_BACKWARD_WISDOM_FAILED -6
#define WISDOM_EXPORT_FILE_OPEN_FAILED -7
#define WISDOM_IMPORT_FILE_OPEN_FAILED -8

#ifndef MIN_SIZE
#define MIN_SIZE 2
#endif

#ifndef MAX_SIZE
#define MAX_SIZE 256
#endif

#ifndef IN_PLACE
#define IN_PLACE true
#endif

#ifndef MAX_ITER
#define MAX_ITER 3
#endif

#ifndef WISDOM_FILE
#define WISDOM_FILE "wisdom.txt"
#endif

int computeEffGrid(int minSize, int maxSize);

#endif /* RANKFFTW_H_ */
