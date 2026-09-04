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

#pragma once

#include <array>

namespace f2dock {

// Configuration for the post-FFT ICP refinement pass (Phase 4 task 4).
//
// The pass runs point-to-plane ICP for each of the top-N final poses,
// using the receptor as the static model cloud and the ligand transformed
// by the candidate 6-DoF pose as the moving cloud. The pose is updated in
// place by composing the ICP delta on top of the FFT-derived rotation /
// translation.
struct IcpRefineConfig {
  bool enabled = false;
  int top_n = 10;
  double inlier_distance = 5.0; // angstroms; 0 disables outlier rejection
  int num_neighbors = 10;       // model PCA neighborhood for normals
};

struct IcpRefineResult {
  bool refined = false;
  double rmsd_before = 0.0;
  double rmsd_after = 0.0;
};

// Refine a single 6-DoF pose against a static receptor point cloud.
//
// `rotation` is row-major: rotation[i][j] is the (i,j) entry of R. The pose
// maps original ligand coordinates into world space via
//     world_ligand = R * orig_ligand + translation.
// On success the rotation/translation are overwritten with the composed pose
//     pose_new = icp_delta * pose_old
// and `refined` is set on the returned result. On failure (cloud too small,
// ICP unable to converge) the inputs are left untouched and `refined` is
// false.
IcpRefineResult refine_pose_point_to_plane(
    const double *receptor_x, const double *receptor_y,
    const double *receptor_z, int n_receptor, const double *ligand_x,
    const double *ligand_y, const double *ligand_z, int n_ligand,
    std::array<std::array<double, 3>, 3> *rotation,
    std::array<double, 3> *translation, const IcpRefineConfig &config);

} // namespace f2dock
