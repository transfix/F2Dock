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

// FFTW backend -- OPT-IN ONLY.
//
// FFTW is GPL-2.0-or-later. Building this file links it, and the resulting
// F2Dock binary is then a combined work covered by the GPL, NOT by F2Dock's
// own LGPL-2.1. It is compiled only when the build is configured with
// -DF2DOCK_FFT_PROVIDER=fftw, which also prints a warning at configure time
// and makes the binary say so at startup.
//
// It is kept because FFTW is the fastest option here (1.41-1.73x over
// PocketFFT on throughput at equal core count) and because reproducing an
// older result may require the exact transform that produced it.

#include "fft-utils/fft_provider.h"

#include <fftw3.h>

#ifdef FFTW_SINGLE_PRECISION
#define FFTW_(sym) fftwf_##sym
#else
#define FFTW_(sym) fftw_##sym
#endif

namespace cvc {
namespace fft {

// FFTW already has a plan object; this is a straight pass-through, so the
// provider seam costs nothing when FFTW is the one selected.
struct plan_impl {
  FFTW_(plan) p;
};

void *alloc(std::size_t bytes) { return FFTW_(malloc)(bytes); }
void release(void *p) { FFTW_(free)(p); }

static plan wrap(FFTW_(plan) p) {
  plan_impl *w = new plan_impl;
  w->p = p;
  return w;
}

plan dft_1d(int n, complex *in, complex *out, int sign, unsigned flags) {
  return wrap(FFTW_(plan_dft_1d)(n, reinterpret_cast<FFTW_(complex) *>(in),
                                 reinterpret_cast<FFTW_(complex) *>(out), sign,
                                 flags));
}

plan dft_2d(int n0, int n1, complex *in, complex *out, int sign,
            unsigned flags) {
  return wrap(FFTW_(plan_dft_2d)(n0, n1, reinterpret_cast<FFTW_(complex) *>(in),
                                 reinterpret_cast<FFTW_(complex) *>(out), sign,
                                 flags));
}

plan dft_3d(int n0, int n1, int n2, complex *in, complex *out, int sign,
            unsigned flags) {
  return wrap(
      FFTW_(plan_dft_3d)(n0, n1, n2, reinterpret_cast<FFTW_(complex) *>(in),
                         reinterpret_cast<FFTW_(complex) *>(out), sign, flags));
}

plan dft_r2c_3d(int n0, int n1, int n2, real *in, complex *out,
                unsigned flags) {
  return wrap(FFTW_(plan_dft_r2c_3d)(
      n0, n1, n2, in, reinterpret_cast<FFTW_(complex) *>(out), flags));
}

plan dft_c2r_3d(int n0, int n1, int n2, complex *in, real *out,
                unsigned flags) {
  return wrap(FFTW_(plan_dft_c2r_3d)(
      n0, n1, n2, reinterpret_cast<FFTW_(complex) *>(in), out, flags));
}

void execute(plan p) {
  if (p)
    FFTW_(execute)(p->p);
}

void destroy(plan p) {
  if (!p)
    return;
  FFTW_(destroy_plan)(p->p);
  delete p;
}

const char *provider_name() { return "FFTW3 (GPL-2.0-or-later)"; }

bool provider_is_gpl() { return true; }

} // namespace fft
} // namespace cvc
