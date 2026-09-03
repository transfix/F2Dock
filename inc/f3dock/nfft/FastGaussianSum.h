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

// FastGaussianSum.h
//
// Phase 5 task 2: a thin wrapper that evaluates a 3D Gaussian sum
//
//     f(y_j) = sum_{k=0..N-1} alpha_k * exp(-c * |y_j - x_k|^2)
//
// at a batch of target points {y_j}, using either an NFFT3 fastsum-backed
// O((N + M) log(N + M)) fast path (when the build was configured with
// `F2DOCK_HAVE_NFFT`) or a deterministic direct O(N*M) reference.
//
// This is the building block that subsequent PRs will route the F2Dock
// electrostatics gridding through under the `nfftAccelerate` flag. The
// direct evaluator is always available and is used as ground truth in the
// unit tests; the NFFT-backed evaluator is wired in a follow-up patch
// once a CI environment with libnfft3 is available.

#ifndef F3DOCK_NFFT_FAST_GAUSSIAN_SUM_H
#define F3DOCK_NFFT_FAST_GAUSSIAN_SUM_H

#include <cstddef>

namespace f3dock {
namespace nfft {

class FastGaussianSum {
public:
  // Construct a Gaussian-sum evaluator.
  //
  //   c           : kernel exponent coefficient; the kernel is
  //                 K(r) = exp(-c * r^2). Must be strictly positive.
  //   num_sources : number of source points (N).
  //   sources     : pointer to N * 3 doubles, interleaved as x0,y0,z0,
  //                 x1,y1,z1, ... Must be non-null when num_sources > 0.
  //   alphas      : pointer to N doubles, the source weights alpha_k.
  //                 Must be non-null when num_sources > 0.
  //   num_targets : number of target points (M).
  //   targets     : pointer to M * 3 doubles, interleaved as above. Must
  //                 be non-null when num_targets > 0.
  //
  // The constructor copies the input pointers; the caller retains ownership
  // and must keep the arrays alive for the lifetime of the object.
  FastGaussianSum(double c, std::size_t num_sources, const double *sources,
                  const double *alphas, std::size_t num_targets,
                  const double *targets);

  ~FastGaussianSum();

  FastGaussianSum(const FastGaussianSum &) = delete;
  FastGaussianSum &operator=(const FastGaussianSum &) = delete;

  // Evaluate f(y_j) for all target points j in [0, num_targets) and write
  // the result into out[0..num_targets). `out` must be non-null and point
  // to at least num_targets doubles.
  //
  // Uses the NFFT3 fastsum fast path when `is_accelerated()` is true and
  // delegates to `evaluate_direct` otherwise.
  void evaluate(double *out) const;

  // Direct O(N*M) reference evaluator. Always available, used as ground
  // truth in unit tests.
  void evaluate_direct(double *out) const;

  // True iff this build was configured with NFFT3 support AND the NFFT
  // fast path has been wired into `evaluate`. Currently always false; the
  // hookup lands in a follow-up patch.
  static bool is_accelerated() noexcept;

  // Inspectors.
  double c() const noexcept { return c_; }
  std::size_t num_sources() const noexcept { return num_sources_; }
  std::size_t num_targets() const noexcept { return num_targets_; }

private:
  double c_;
  std::size_t num_sources_;
  const double *sources_;
  const double *alphas_;
  std::size_t num_targets_;
  const double *targets_;
};

} // namespace nfft
} // namespace f3dock

#endif // F3DOCK_NFFT_FAST_GAUSSIAN_SUM_H
