#pragma once

#include "f3dock/flex/FlexKinematics.h"
#include "f3dock/flex/FlexSampler.h"

#include <cstddef>

namespace f3dock {
namespace flex {

// Apply every hinge then every shear in `state` to the flat coordinate
// arrays (x, y, z), each containing `n` doubles. Operates in place.
//
// F2Dock stores receptor and ligand atom coordinates as three separate
// double* arrays (xkAOrig, ykAOrig, zkAOrig). This helper bridges that
// layout to the std::vector<Point3>-based FlexKinematics API without
// forcing callers to repack manually.
//
// When `state` has no hinges and no shears (the identity state at
// index 0 of an unconfigured FlexSamplingConfig) this is a no-op and
// the arrays are left untouched.
//
// Returns true on success. Returns false (leaving the arrays unchanged
// from before the failing transform) if:
//   * any of x, y, z is nullptr while n > 0,
//   * any hinge has a zero-length axis_direction,
//   * any shear has a zero-length plane_normal.
bool apply_flex_state(const FlexState &state, double *x, double *y, double *z,
                      std::size_t n);

} // namespace flex
} // namespace f3dock
