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

/*
  Link-time stand-in for the GPL sparse 3-D FFT module, compiled in place of
  sparsefft3.cpp / sparsefft3-plan.cpp when F2DOCK_ENABLE_GPL_COMPONENTS is
  OFF (the default). See inc/fft-utils/sparsefft3-api.h.

  Nothing here is reachable in a well-formed run: the engine forces
  useSparseFFT off when the module is absent, which pins both plan pointers
  to NULL and routes every transform through dense FFTW. The bodies exist so
  that a future caller which forgets that guard fails loudly rather than
  silently producing a wrong score.
*/

#include "fft-utils/sparsefft3-api.h"

#ifndef F2DOCK_WITH_SPARSEFFT

#include <stdio.h>
#include <stdlib.h>

static void sparsefft3_unavailable(const char *fn) {
  fprintf(stderr,
          "F2Dock: %s was called, but the sparse 3-D FFT module is not part "
          "of this build.\n"
          "        It is GPL-licensed and is only compiled when F2Dock is "
          "configured with\n"
          "        -DF2DOCK_ENABLE_GPL_COMPONENTS=ON. Set 'useSparseFFT "
          "false' in the\n"
          "        parameter file to use the dense FFTW path instead.\n",
          fn);
}

sparse3DFFT_plan
sparse3DFFT_create_plan(int nx, int ny, int nz, int dir, int flags,
                        sparse3DFFT_sparsedir sparsedir,
                        sparse3DFFT_nonzero_func nonzero, void *nonzero_data,
                        FFTW_complex *data_in, FFTW_complex *data_out) {
  (void)nx;
  (void)ny;
  (void)nz;
  (void)dir;
  (void)flags;
  (void)sparsedir;
  (void)nonzero;
  (void)nonzero_data;
  (void)data_in;
  (void)data_out;
  sparsefft3_unavailable("sparse3DFFT_create_plan");
  return NULL;
}

void sparse3DFFT_destroy_plan(sparse3DFFT_plan p) { (void)p; }

void sparse3DFFT(sparse3DFFT_plan p, FFTW_complex *data_in,
                 FFTW_complex *data_out) {
  (void)p;
  (void)data_in;
  (void)data_out;
  sparsefft3_unavailable("sparse3DFFT");
  abort();
}

#endif /* !F2DOCK_WITH_SPARSEFFT */
