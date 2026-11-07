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

#ifndef F2DOCK_FFT_UTILS_FFT_PROVIDER_H
#define F2DOCK_FFT_UTILS_FFT_PROVIDER_H

/*
  Provider-neutral FFT interface.

  WHY THIS EXISTS. FFTW is GPL-2.0-or-later. F2Dock is LGPL-2.1-only, so a
  binary linking FFTW is a combined work under the GPL -- which is the whole
  reason this seam exists, not performance. Gating the vendored GPL modules
  (sparsefft3, libicp) cleaned the source tree but left the linked library
  problem untouched. See THIRD-PARTY.md.

  WHY IT IS SMALL. The engine never touched FFTW directly: every call already
  went through the FFTW_* macros in fftwPrecision.h. The surface those macros
  cover, in the core engine, is only what appears below. `plan_many_dft` and
  `execute_dft` are used solely by the (removed) sparsefft3 module, and FFTW's
  wisdom API solely by rank-fftw's grid-timing utility, so neither is part of
  this interface.

  SEMANTICS. Deliberately FFTW's, so a provider swap is not a behaviour change:

    * Transforms are UNNORMALISED in both directions, exactly like FFTW. A
      forward followed by a backward scales by n0*n1*n2. Docking.cpp relies on
      this; a provider that silently normalised would change every score.
    * sign is FFTW_FORWARD (-1) or FFTW_BACKWARD (+1).
    * r2c/c2r use FFTW's half-complex layout: the complex side has
      n0 * n1 * (n2/2 + 1) elements.
    * A plan captures its buffers at creation. execute() re-runs on those same
      buffers, which is how the engine reuses one plan across many rotations.
    * flags are advisory. FFTW uses them to pick an algorithm; providers with
      no planning stage ignore them.

  ONE DELIBERATE DIVERGENCE: FFTW's c2r destroys its input array. Providers
  here are not required to, and PocketFFT does not. Nothing in F2Dock reads a
  c2r input after transforming it, so preserving it is strictly safer.
*/

#include <cstddef>

// Precision follows the same switch the rest of the tree uses.
#ifdef FFTW_SINGLE_PRECISION
#define F2DOCK_FFT_REAL float
#else
#define F2DOCK_FFT_REAL double
#endif

namespace cvc {
namespace fft {

typedef F2DOCK_FFT_REAL real;
// Interleaved re/im, layout-compatible with fftw_complex and
// std::complex<real>, so buffers pass between providers without a copy.
typedef real complex[2];

// Opaque; each provider defines its own.
struct plan_impl;
typedef plan_impl *plan;

// Sign and flag values are FFTW's, so a build configured against FFTW can map
// straight through without translating them.
enum {
  sign_forward = -1,
  sign_backward = 1,
};

// Aligned allocation suitable for the provider's vectorised kernels.
void *alloc(std::size_t bytes);
void release(void *p);

plan dft_1d(int n, complex *in, complex *out, int sign, unsigned flags);
plan dft_2d(int n0, int n1, complex *in, complex *out, int sign,
            unsigned flags);
plan dft_3d(int n0, int n1, int n2, complex *in, complex *out, int sign,
            unsigned flags);
plan dft_r2c_3d(int n0, int n1, int n2, real *in, complex *out, unsigned flags);
plan dft_c2r_3d(int n0, int n1, int n2, complex *in, real *out, unsigned flags);

void execute(plan p);
void destroy(plan p);

// Human-readable provider name, for the startup banner so a run log records
// which transform produced its numbers.
const char *provider_name();

// True when the linked provider is GPL (i.e. FFTW), so callers can warn.
bool provider_is_gpl();

} // namespace fft
} // namespace cvc

#endif // F2DOCK_FFT_UTILS_FFT_PROVIDER_H
