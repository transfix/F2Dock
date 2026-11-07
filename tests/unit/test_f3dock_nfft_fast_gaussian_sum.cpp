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

// Unit tests for f3dock::nfft::FastGaussianSum.
//
// Locks in:
//   * direct evaluator matches a hand-computed reference,
//   * `evaluate()` matches `evaluate_direct()` exactly when the build is
//     not NFFT-accelerated (which is the case in every current CI runner),
//   * input validation rejects bad arguments,
//   * single-point identity: f(x_k) >= alpha_k (when alphas are positive).

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "f3dock/nfft/FastGaussianSum.h"

namespace {

// Compute a single Gaussian sum value via a hand-written reference loop.
double reference_sum(double c, const std::vector<double> &sources,
                     const std::vector<double> &alphas,
                     const std::array<double, 3> &target) {
  const std::size_t n = alphas.size();
  double acc = 0.0;
  for (std::size_t k = 0; k < n; ++k) {
    const double dx = target[0] - sources[3 * k + 0];
    const double dy = target[1] - sources[3 * k + 1];
    const double dz = target[2] - sources[3 * k + 2];
    const double r2 = dx * dx + dy * dy + dz * dz;
    acc += alphas[k] * std::exp(-c * r2);
  }
  return acc;
}

} // namespace

TEST(FastGaussianSum, ConstructorRejectsNonPositiveC) {
  std::vector<double> srcs = {0.0, 0.0, 0.0};
  std::vector<double> alphas = {1.0};
  std::vector<double> tgts = {1.0, 0.0, 0.0};
  EXPECT_THROW(f3dock::nfft::FastGaussianSum(0.0, 1, srcs.data(), alphas.data(),
                                             1, tgts.data()),
               std::invalid_argument);
  EXPECT_THROW(f3dock::nfft::FastGaussianSum(-1.0, 1, srcs.data(),
                                             alphas.data(), 1, tgts.data()),
               std::invalid_argument);
}

TEST(FastGaussianSum, ConstructorRejectsNullSourcesWhenNPositive) {
  std::vector<double> alphas = {1.0};
  std::vector<double> tgts = {1.0, 0.0, 0.0};
  EXPECT_THROW(f3dock::nfft::FastGaussianSum(1.0, 1, nullptr, alphas.data(), 1,
                                             tgts.data()),
               std::invalid_argument);
}

TEST(FastGaussianSum, EvaluateDirectMatchesReference) {
  const double c = 2.5;
  std::vector<double> srcs = {
      0.0, 0.0, 0.0, 1.0, 0.5, -0.25, -0.5, 0.75, 0.5,
  };
  std::vector<double> alphas = {1.0, -2.0, 0.7};
  std::vector<double> tgts = {
      0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 2.0, -1.0, 0.5,
  };
  const std::size_t M = tgts.size() / 3;
  f3dock::nfft::FastGaussianSum sum(c, alphas.size(), srcs.data(),
                                    alphas.data(), M, tgts.data());
  std::vector<double> out(M, std::nan(""));
  sum.evaluate_direct(out.data());
  for (std::size_t j = 0; j < M; ++j) {
    const std::array<double, 3> y{tgts[3 * j + 0], tgts[3 * j + 1],
                                  tgts[3 * j + 2]};
    EXPECT_NEAR(out[j], reference_sum(c, srcs, alphas, y), 1e-12) << "j=" << j;
  }
}

TEST(FastGaussianSum, EvaluateMatchesDirectWhenUnaccelerated) {
  const double c = 1.0;
  std::vector<double> srcs = {0.1, -0.2, 0.3, 0.4, 0.5, -0.6};
  std::vector<double> alphas = {0.5, 1.5};
  std::vector<double> tgts = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, -0.5, 0.25, 0.0};
  const std::size_t M = tgts.size() / 3;
  f3dock::nfft::FastGaussianSum sum(c, alphas.size(), srcs.data(),
                                    alphas.data(), M, tgts.data());
  std::vector<double> direct(M), fast(M);
  sum.evaluate_direct(direct.data());
  sum.evaluate(fast.data());
  // When `is_accelerated()` is false, `evaluate` must be bit-identical to
  // `evaluate_direct` (it's the same routine called twice).
  if (!f3dock::nfft::FastGaussianSum::is_accelerated()) {
    for (std::size_t j = 0; j < M; ++j) {
      EXPECT_EQ(fast[j], direct[j]) << "j=" << j;
    }
  } else {
    // When the NFFT3 fast path is eventually wired up, allow a small
    // discrepancy bounded by the fastsum regularization tolerance.
    for (std::size_t j = 0; j < M; ++j) {
      EXPECT_NEAR(fast[j], direct[j], 1e-6) << "j=" << j;
    }
  }
}

TEST(FastGaussianSum, SinglePointIdentityIsBoundedBelow) {
  // f(x_0) = alpha_0 * exp(0) + sum_{k!=0} alpha_k * exp(-c * |x_0 - x_k|^2)
  // With positive alphas, f(x_0) >= alpha_0.
  const double c = 4.0;
  std::vector<double> srcs = {0.0, 0.0, 0.0, 2.0, 2.0, 2.0};
  std::vector<double> alphas = {3.0, 1.0};
  std::vector<double> tgts = {0.0, 0.0, 0.0};
  f3dock::nfft::FastGaussianSum sum(c, 2, srcs.data(), alphas.data(), 1,
                                    tgts.data());
  std::vector<double> out(1);
  sum.evaluate_direct(out.data());
  EXPECT_GE(out[0], 3.0);
  // Upper bound: alpha_0 * 1 + alpha_1 * 1 = 4
  EXPECT_LE(out[0], 4.0);
}

TEST(FastGaussianSum, EmptySourcesProducesZero) {
  std::vector<double> tgts = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  f3dock::nfft::FastGaussianSum sum(1.0, 0, nullptr, nullptr, 2, tgts.data());
  std::vector<double> out = {7.0, 9.0};
  sum.evaluate_direct(out.data());
  EXPECT_EQ(out[0], 0.0);
  EXPECT_EQ(out[1], 0.0);
}

TEST(FastGaussianSum, Inspectors) {
  std::vector<double> srcs = {0.0, 0.0, 0.0};
  std::vector<double> alphas = {1.0};
  std::vector<double> tgts = {1.0, 0.0, 0.0};
  f3dock::nfft::FastGaussianSum sum(2.0, 1, srcs.data(), alphas.data(), 1,
                                    tgts.data());
  EXPECT_DOUBLE_EQ(sum.c(), 2.0);
  EXPECT_EQ(sum.num_sources(), 1u);
  EXPECT_EQ(sum.num_targets(), 1u);
}
