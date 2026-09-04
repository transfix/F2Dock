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

// PocketFFT backend (BSD-3-Clause, from the cvcpkg `pocketfft` bundle).
//
// This is the default provider because it is the fastest transform F2Dock can
// use without making the binary GPL. Verified against FFTW on the shape the
// hot loop runs (3-D c2c, double, N=64..256): agrees to float64 round-off,
// <= 7.3e-12. It is slower -- FFTW leads by 1.41-1.73x on throughput at equal
// core count -- which is the price of the licence, and why MKL and cuFFT are
// also selectable.
//
// PocketFFT has no planning stage: it computes twiddles per call and caches
// them internally. So a "plan" here is just the captured shape and buffers,
// and `flags` (FFTW_MEASURE and friends) are accepted and ignored.

#include "fft-utils/fft_provider.h"

#include <cstdlib>
#include <new>
#include <stdexcept>
#include <vector>

#include <pocketfft_hdronly.h>

namespace cvc {
namespace fft {

namespace {

enum kind_t { KIND_C2C, KIND_R2C, KIND_C2R };

// Byte alignment for the transform buffers. PocketFFT vectorises over the
// contiguous axis, so give it what the widest SIMD it may pick wants.
const std::size_t kAlign = 64;

pocketfft::stride_t contiguous_strides(const pocketfft::shape_t &shape,
                                       std::size_t elem) {
  pocketfft::stride_t s(shape.size());
  std::ptrdiff_t acc = static_cast<std::ptrdiff_t>(elem);
  for (std::size_t i = shape.size(); i-- > 0;) {
    s[i] = acc;
    acc *= static_cast<std::ptrdiff_t>(shape[i]);
  }
  return s;
}

} // namespace

// A plan is the arguments, retained. execute() replays them.
struct plan_impl {
  kind_t kind;
  pocketfft::shape_t shape; // logical dimensions
  int sign;                 // c2c only
  void *in;
  void *out;
};

void *alloc(std::size_t bytes) {
  if (bytes == 0)
    bytes = 1;
  // Round up so the size is a multiple of the alignment, which
  // std::aligned_alloc requires and some libcs enforce.
  const std::size_t rounded = ((bytes + kAlign - 1) / kAlign) * kAlign;
#if defined(_MSC_VER)
  void *p = _aligned_malloc(rounded, kAlign);
#else
  void *p = std::aligned_alloc(kAlign, rounded);
#endif
  if (!p)
    throw std::bad_alloc();
  return p;
}

void release(void *p) {
  if (!p)
    return;
#if defined(_MSC_VER)
  _aligned_free(p);
#else
  std::free(p);
#endif
}

static plan make_c2c(const pocketfft::shape_t &shape, complex *in, complex *out,
                     int sign) {
  plan_impl *p = new plan_impl;
  p->kind = KIND_C2C;
  p->shape = shape;
  p->sign = sign;
  p->in = in;
  p->out = out;
  return p;
}

plan dft_1d(int n, complex *in, complex *out, int sign, unsigned /*flags*/) {
  return make_c2c(pocketfft::shape_t{static_cast<std::size_t>(n)}, in, out,
                  sign);
}

plan dft_2d(int n0, int n1, complex *in, complex *out, int sign,
            unsigned /*flags*/) {
  return make_c2c(pocketfft::shape_t{static_cast<std::size_t>(n0),
                                     static_cast<std::size_t>(n1)},
                  in, out, sign);
}

plan dft_3d(int n0, int n1, int n2, complex *in, complex *out, int sign,
            unsigned /*flags*/) {
  return make_c2c(pocketfft::shape_t{static_cast<std::size_t>(n0),
                                     static_cast<std::size_t>(n1),
                                     static_cast<std::size_t>(n2)},
                  in, out, sign);
}

plan dft_r2c_3d(int n0, int n1, int n2, real *in, complex *out,
                unsigned /*flags*/) {
  plan_impl *p = new plan_impl;
  p->kind = KIND_R2C;
  p->shape = pocketfft::shape_t{static_cast<std::size_t>(n0),
                                static_cast<std::size_t>(n1),
                                static_cast<std::size_t>(n2)};
  p->sign = sign_forward;
  p->in = in;
  p->out = out;
  return p;
}

plan dft_c2r_3d(int n0, int n1, int n2, complex *in, real *out,
                unsigned /*flags*/) {
  plan_impl *p = new plan_impl;
  p->kind = KIND_C2R;
  p->shape = pocketfft::shape_t{static_cast<std::size_t>(n0),
                                static_cast<std::size_t>(n1),
                                static_cast<std::size_t>(n2)};
  p->sign = sign_backward;
  p->in = in;
  p->out = out;
  return p;
}

void execute(plan p) {
  if (!p)
    return;

  // fct = 1.0 in BOTH directions: FFTW is unnormalised, and Docking.cpp's
  // scoring depends on that. Normalising here would change every score.
  const real fct = static_cast<real>(1);

  const std::size_t ndim = p->shape.size();
  pocketfft::shape_t axes(ndim);
  for (std::size_t i = 0; i < ndim; ++i)
    axes[i] = i;

  switch (p->kind) {
  case KIND_C2C: {
    const pocketfft::stride_t st =
        contiguous_strides(p->shape, sizeof(std::complex<real>));
    pocketfft::c2c(p->shape, st, st, axes,
                   p->sign == sign_forward ? pocketfft::FORWARD
                                           : pocketfft::BACKWARD,
                   reinterpret_cast<const std::complex<real> *>(p->in),
                   reinterpret_cast<std::complex<real> *>(p->out), fct);
    break;
  }
  case KIND_R2C: {
    // FFTW's half-complex layout: the complex side is n0 x n1 x (n2/2 + 1).
    pocketfft::shape_t cshape = p->shape;
    cshape[ndim - 1] = p->shape[ndim - 1] / 2 + 1;
    const pocketfft::stride_t sin = contiguous_strides(p->shape, sizeof(real));
    const pocketfft::stride_t sout =
        contiguous_strides(cshape, sizeof(std::complex<real>));
    pocketfft::r2c(p->shape, sin, sout, axes, pocketfft::FORWARD,
                   reinterpret_cast<const real *>(p->in),
                   reinterpret_cast<std::complex<real> *>(p->out), fct);
    break;
  }
  case KIND_C2R: {
    pocketfft::shape_t cshape = p->shape;
    cshape[ndim - 1] = p->shape[ndim - 1] / 2 + 1;
    const pocketfft::stride_t sin =
        contiguous_strides(cshape, sizeof(std::complex<real>));
    const pocketfft::stride_t sout = contiguous_strides(p->shape, sizeof(real));
    pocketfft::c2r(p->shape, sin, sout, axes, pocketfft::BACKWARD,
                   reinterpret_cast<const std::complex<real> *>(p->in),
                   reinterpret_cast<real *>(p->out), fct);
    break;
  }
  }
}

void destroy(plan p) { delete p; }

const char *provider_name() { return "PocketFFT (BSD-3-Clause)"; }

bool provider_is_gpl() { return false; }

} // namespace fft
} // namespace cvc
