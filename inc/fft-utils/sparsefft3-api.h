/*
  Copyright 2026 The University of Texas at Austin

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

#ifndef F2DOCK_FFT_UTILS_SPARSEFFT3_API_H
#define F2DOCK_FFT_UTILS_SPARSEFFT3_API_H

/*
  Licence seam for the optional sparse 3-D FFT module.

  sparsefft3.h / sparsefft3.cpp / sparsefft3-plan.cpp are GPL-2.0-or-later
  (MIT 2000, modified at CVC Lab), while F2Dock as a whole is LGPL-2.1-only.
  The LGPL core therefore must not #include the GPL headers directly: they
  carry substantive definitions -- structs, enums, macros -- not just
  prototypes.

  This header is the seam. Include it instead of sparsefft3.h. When the build
  opts in to the GPL components (F2DOCK_ENABLE_GPL_COMPONENTS=ON, which sets
  F2DOCK_WITH_SPARSEFFT), it forwards to the real header. Otherwise it
  declares the same entry points against an incomplete plan type, which is
  all the LGPL core ever needs: the plan is only ever held and passed as an
  opaque pointer.

  In the default build the plan pointers are held at NULL (Docking.cpp, at
  the `if (!useSparseFFT)` branch) and every call site falls through to the
  dense FFTW path, so none of the functions declared below is reached.
  src/fft-utils/sparsefft3-stub.cpp satisfies the link.
*/

#include "fftw3.h"
#include "fftwPrecision.h"

#ifdef F2DOCK_WITH_SPARSEFFT

#include "sparsefft3.h"

#else

/* Incomplete type: never dereferenced outside the GPL module. */
typedef struct sparse3DFFT_plan_opaque_ *sparse3DFFT_plan;

typedef enum {
  SPARSE3DFFT_SPARSEINPUT,
  SPARSE3DFFT_SPARSEOUTPUT
} sparse3DFFT_sparsedir;

typedef int (*sparse3DFFT_nonzero_func)(int x[3], void *data);

extern sparse3DFFT_plan
sparse3DFFT_create_plan(int nx, int ny, int nz, int dir, int flags,
                        sparse3DFFT_sparsedir sparsedir,
                        sparse3DFFT_nonzero_func nonzero, void *nonzero_data,
                        FFTW_complex *data_in, FFTW_complex *data_out);

extern void sparse3DFFT_destroy_plan(sparse3DFFT_plan p);

extern void sparse3DFFT(sparse3DFFT_plan p, FFTW_complex *data_in,
                        FFTW_complex *data_out);

#endif /* F2DOCK_WITH_SPARSEFFT */

#endif /* F2DOCK_FFT_UTILS_SPARSEFFT3_API_H */
