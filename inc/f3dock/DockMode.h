// SPDX-License-Identifier: LGPL-2.1-only
//
// DockMode.h
//
// Runtime selector between the two docking pipelines that share the
// consolidated F2Dock executable:
//
//   F2 mode (kF2Dock) -- "Fast Fourier" rigid-body docking. The
//                        original F2Dock pipeline. Both partners are
//                        treated as rigid bodies and the FFT-accelerated
//                        scoring is applied directly to the input
//                        conformations.
//
//   F3 mode (kF3Dock) -- "Flexible Fast Fourier" docking. Adds the
//                        F3Dock-derived sampling stages (tripeptide loop
//                        closure, hinge / shear protein flexibility,
//                        domain-graph composition, ICP refinement) on
//                        top of the F2Dock FFT scoring core.
//
// The mode is a runtime decision so a single binary can serve both
// pipelines. F3-specific modules under inc/f3dock/{loop_closure,flex,
// icp,domain} are only invoked when the active mode is kF3Dock.
//
// Selecting the mode:
//
//   * CLI:    F2Dock --mode=f2|f3  parameterFile
//             F2Dock --f2dock      parameterFile
//             F2Dock --f3dock      parameterFile
//
//   * Input:  add 'dockMode f2' or 'dockMode f3' to the .inp file.
//             Synonyms accepted: 'rigid' / 'flexible'.
//
//   * Default: kF2Dock (rigid). Existing F2Dock workflows keep working
//              unchanged with no input changes.
//
// CLI overrides input-file value when both are present.

#ifndef F3DOCK_DOCKMODE_H
#define F3DOCK_DOCKMODE_H

#include <cstring>

namespace f3dock {

enum class DockMode {
  kF2Dock = 0, // rigid Fast Fourier docking (legacy F2Dock)
  kF3Dock = 1  // Flexible Fast Fourier docking (F2Dock + F3Dock stages)
};

inline const char *to_string(DockMode m) {
  return (m == DockMode::kF3Dock) ? "F3Dock (flexible)" : "F2Dock (rigid)";
}

// Parse a mode token from a CLI flag value or input-file value.
// Accepted (case-insensitive): "f2", "f2dock", "rigid"
//                              "f3", "f3dock", "flex", "flexible"
// Returns true on success and writes to *out; false on unrecognized input.
inline bool parse_dock_mode(const char *s, DockMode *out) {
  if (s == nullptr || out == nullptr)
    return false;
  if (strcasecmp(s, "f2") == 0 || strcasecmp(s, "f2dock") == 0 ||
      strcasecmp(s, "rigid") == 0) {
    *out = DockMode::kF2Dock;
    return true;
  }
  if (strcasecmp(s, "f3") == 0 || strcasecmp(s, "f3dock") == 0 ||
      strcasecmp(s, "flex") == 0 || strcasecmp(s, "flexible") == 0) {
    *out = DockMode::kF3Dock;
    return true;
  }
  return false;
}

} // namespace f3dock

#endif // F3DOCK_DOCKMODE_H
