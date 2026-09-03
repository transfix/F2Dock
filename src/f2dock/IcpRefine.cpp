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

#include "f2dock/IcpRefine.h"

#include <algorithm>
#include <vector>

#include "f3dock/icp/IcpAligner.h"

namespace f2dock {

using f3dock::icp::IcpAligner;
using f3dock::icp::IcpAlignmentResult;
using f3dock::icp::Point3;

IcpRefineResult refine_pose_point_to_plane(
    const double *receptor_x, const double *receptor_y,
    const double *receptor_z, int n_receptor, const double *ligand_x,
    const double *ligand_y, const double *ligand_z, int n_ligand,
    std::array<std::array<double, 3>, 3> *rotation,
    std::array<double, 3> *translation, const IcpRefineConfig &config) {
  IcpRefineResult out;

  if (!rotation || !translation || !receptor_x || !receptor_y || !receptor_z ||
      !ligand_x || !ligand_y || !ligand_z) {
    return out;
  }
  if (n_receptor < IcpAligner::kMinPoints ||
      n_ligand < IcpAligner::kMinPoints) {
    return out;
  }

  std::vector<Point3> model;
  model.reserve(static_cast<size_t>(n_receptor));
  for (int i = 0; i < n_receptor; ++i) {
    model.push_back({receptor_x[i], receptor_y[i], receptor_z[i]});
  }

  const auto &R = *rotation;
  const auto &t = *translation;

  std::vector<Point3> moving;
  moving.reserve(static_cast<size_t>(n_ligand));
  for (int i = 0; i < n_ligand; ++i) {
    const double x = ligand_x[i];
    const double y = ligand_y[i];
    const double z = ligand_z[i];
    moving.push_back({
        R[0][0] * x + R[0][1] * y + R[0][2] * z + t[0],
        R[1][0] * x + R[1][1] * y + R[1][2] * z + t[1],
        R[2][0] * x + R[2][1] * y + R[2][2] * z + t[2],
    });
  }

  int nn = config.num_neighbors;
  if (nn < 3) {
    nn = 3;
  }
  if (static_cast<size_t>(nn) > model.size()) {
    nn = static_cast<int>(model.size());
  }

  IcpAligner aligner;
  std::vector<Point3> aligned;
  IcpAlignmentResult result;
  const bool ok = aligner.alignPointToPlane(model, moving, aligned, &result,
                                            config.inlier_distance, nn);
  if (!ok) {
    return out;
  }

  // Compose the ICP delta on top of the existing pose:
  //   pose_new = icp_delta * pose_old
  // i.e. R_new = R_d * R, t_new = R_d * t + t_d.
  const auto &Rd = result.rotation;
  const auto &td = result.translation;

  std::array<std::array<double, 3>, 3> Rnew{};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      Rnew[i][j] = Rd[i][0] * R[0][j] + Rd[i][1] * R[1][j] + Rd[i][2] * R[2][j];
    }
  }
  std::array<double, 3> tnew{
      Rd[0][0] * t[0] + Rd[0][1] * t[1] + Rd[0][2] * t[2] + td[0],
      Rd[1][0] * t[0] + Rd[1][1] * t[1] + Rd[1][2] * t[2] + td[1],
      Rd[2][0] * t[0] + Rd[2][1] * t[1] + Rd[2][2] * t[2] + td[2],
  };

  *rotation = Rnew;
  *translation = tnew;
  out.refined = true;
  out.rmsd_before = result.rmsd_before;
  out.rmsd_after = result.rmsd_after;
  return out;
}

} // namespace f2dock
