// FastGaussianSum.cpp -- implementation. See FastGaussianSum.h for the
// contract. The NFFT3 fastsum path is currently a no-op stub gated by
// F2DOCK_HAVE_NFFT_WIRED so that the library compiles uniformly on every
// platform in the CI matrix; the real fastsum hookup lands in a follow-up
// patch once libnfft3 is available in CI.

#include "f3dock/nfft/FastGaussianSum.h"

#include <cmath>
#include <stdexcept>

namespace f3dock {
namespace nfft {

FastGaussianSum::FastGaussianSum(double c, std::size_t num_sources,
                                 const double *sources, const double *alphas,
                                 std::size_t num_targets, const double *targets)
    : c_(c), num_sources_(num_sources), sources_(sources), alphas_(alphas),
      num_targets_(num_targets), targets_(targets) {
  if (!(c_ > 0.0)) {
    throw std::invalid_argument(
        "FastGaussianSum: kernel coefficient c must be > 0");
  }
  if (num_sources_ > 0 && (sources_ == nullptr || alphas_ == nullptr)) {
    throw std::invalid_argument(
        "FastGaussianSum: sources/alphas must be non-null when N > 0");
  }
  if (num_targets_ > 0 && targets_ == nullptr) {
    throw std::invalid_argument(
        "FastGaussianSum: targets must be non-null when M > 0");
  }
}

FastGaussianSum::~FastGaussianSum() = default;

void FastGaussianSum::evaluate_direct(double *out) const {
  if (out == nullptr) {
    throw std::invalid_argument(
        "FastGaussianSum::evaluate_direct: out is null");
  }
  for (std::size_t j = 0; j < num_targets_; ++j) {
    const double yx = targets_[3 * j + 0];
    const double yy = targets_[3 * j + 1];
    const double yz = targets_[3 * j + 2];
    double acc = 0.0;
    for (std::size_t k = 0; k < num_sources_; ++k) {
      const double dx = yx - sources_[3 * k + 0];
      const double dy = yy - sources_[3 * k + 1];
      const double dz = yz - sources_[3 * k + 2];
      const double r2 = dx * dx + dy * dy + dz * dz;
      acc += alphas_[k] * std::exp(-c_ * r2);
    }
    out[j] = acc;
  }
}

void FastGaussianSum::evaluate(double *out) const {
  // Placeholder for the NFFT3 fastsum fast path. When wired, this branch
  // will set up an `fastsum_plan`, populate source/target nodes and
  // weights, call `fastsum_trafo`, and copy `plan.f` into `out`. Until
  // libnfft3 is available in the project's CI matrix, the fast path is
  // intentionally a hard fallback to the direct evaluator so that the
  // observable behavior is identical on every platform.
  evaluate_direct(out);
}

bool FastGaussianSum::is_accelerated() noexcept {
  // Will return `true` once the NFFT3 fastsum hookup lands AND the build
  // was configured with `F2DOCK_HAVE_NFFT`. For now, the fast path is
  // deliberately disabled even on builds where NFFT3 is available so the
  // direct fallback is exercised uniformly.
  return false;
}

} // namespace nfft
} // namespace f3dock
