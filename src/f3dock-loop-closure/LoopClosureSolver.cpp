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

#include "f3dock/loop_closure/LoopClosureSolver.h"

#include <cmath>

#include "tripep_closure.h"

namespace f3dock {
namespace loop_closure {
namespace {

enum AtomIndex {
  kN1 = 0,
  kCA1 = 1,
  kC1 = 2,
  kN2 = 3,
  kCA2 = 4,
  kC2 = 5,
  kN3 = 6,
  kCA3 = 7,
  kC3 = 8,
};

Point3 subtract(const Point3 &a, const Point3 &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

double dot(const Point3 &a, const Point3 &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Point3 cross(const Point3 &a, const Point3 &b) {
  return {
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0],
  };
}

double norm(const Point3 &v) { return std::sqrt(dot(v, v)); }

double bondLength(const Point3 &a, const Point3 &b) {
  return norm(subtract(a, b));
}

double bondAngle(const Point3 &a, const Point3 &b, const Point3 &c) {
  const Point3 u = subtract(a, b);
  const Point3 v = subtract(c, b);
  const double nu = norm(u);
  const double nv = norm(v);
  if (nu == 0.0 || nv == 0.0) {
    return 0.0;
  }
  double cos_theta = dot(u, v) / (nu * nv);
  if (cos_theta > 1.0) {
    cos_theta = 1.0;
  }
  if (cos_theta < -1.0) {
    cos_theta = -1.0;
  }
  return std::acos(cos_theta);
}

double torsion(const Point3 &p1, const Point3 &p2, const Point3 &p3,
               const Point3 &p4) {
  const Point3 b1 = subtract(p2, p1);
  const Point3 b2 = subtract(p3, p2);
  const Point3 b3 = subtract(p4, p3);

  const Point3 n1 = cross(b1, b2);
  const Point3 n2 = cross(b2, b3);
  const double n1n = norm(n1);
  const double n2n = norm(n2);
  const double b2n = norm(b2);
  if (n1n == 0.0 || n2n == 0.0 || b2n == 0.0) {
    return 0.0;
  }

  const Point3 m1 = cross(n1, {b2[0] / b2n, b2[1] / b2n, b2[2] / b2n});
  const double x = dot(n1, n2);
  const double y = dot(m1, n2);
  return std::atan2(y, x);
}

void pointToArray(const Point3 &p, double out[3]) {
  out[0] = p[0];
  out[1] = p[1];
  out[2] = p[2];
}

} // namespace

bool LoopClosureSolver::solve(const LoopClosureInput &input,
                              std::vector<Backbone9> &solutions) const {
  const Backbone9 &bb = input.backbone;

  double b_len[6];
  b_len[0] = bondLength(bb[kCA1], bb[kC1]);
  b_len[1] = bondLength(bb[kC1], bb[kN2]);
  b_len[2] = bondLength(bb[kN2], bb[kCA2]);
  b_len[3] = bondLength(bb[kCA2], bb[kC2]);
  b_len[4] = bondLength(bb[kC2], bb[kN3]);
  b_len[5] = bondLength(bb[kN3], bb[kCA3]);

  double b_ang[7];
  b_ang[0] = bondAngle(bb[kN1], bb[kCA1], bb[kC1]);
  b_ang[1] = bondAngle(bb[kCA1], bb[kC1], bb[kN2]);
  b_ang[2] = bondAngle(bb[kC1], bb[kN2], bb[kCA2]);
  b_ang[3] = bondAngle(bb[kN2], bb[kCA2], bb[kC2]);
  b_ang[4] = bondAngle(bb[kCA2], bb[kC2], bb[kN3]);
  b_ang[5] = bondAngle(bb[kC2], bb[kN3], bb[kCA3]);
  b_ang[6] = bondAngle(bb[kN3], bb[kCA3], bb[kC3]);

  double t_ang[2];
  t_ang[0] = torsion(bb[kCA1], bb[kC1], bb[kN2], bb[kCA2]);
  t_ang[1] = torsion(bb[kCA2], bb[kC2], bb[kN3], bb[kCA3]);

  double r_n1[3], r_a1[3], r_a3[3], r_c3[3];
  pointToArray(input.constraints.n1, r_n1);
  pointToArray(input.constraints.ca1, r_a1);
  pointToArray(input.constraints.ca3, r_a3);
  pointToArray(input.constraints.c3, r_c3);

  double r_soln_n[max_soln][3][3] = {};
  double r_soln_a[max_soln][3][3] = {};
  double r_soln_c[max_soln][3][3] = {};
  int n_soln = 0;

  tripep_closure closure;
  closure.initialize_loop_closure(b_len, b_ang, t_ang);
  closure.solve_3pep_poly(r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c,
                          &n_soln);

  solutions.clear();
  for (int i = 0; i < n_soln; ++i) {
    Backbone9 sol{};
    sol[kN1] = {r_soln_n[i][0][0], r_soln_n[i][0][1], r_soln_n[i][0][2]};
    sol[kCA1] = {r_soln_a[i][0][0], r_soln_a[i][0][1], r_soln_a[i][0][2]};
    sol[kC1] = {r_soln_c[i][0][0], r_soln_c[i][0][1], r_soln_c[i][0][2]};
    sol[kN2] = {r_soln_n[i][1][0], r_soln_n[i][1][1], r_soln_n[i][1][2]};
    sol[kCA2] = {r_soln_a[i][1][0], r_soln_a[i][1][1], r_soln_a[i][1][2]};
    sol[kC2] = {r_soln_c[i][1][0], r_soln_c[i][1][1], r_soln_c[i][1][2]};
    sol[kN3] = {r_soln_n[i][2][0], r_soln_n[i][2][1], r_soln_n[i][2][2]};
    sol[kCA3] = {r_soln_a[i][2][0], r_soln_a[i][2][1], r_soln_a[i][2][2]};
    sol[kC3] = {r_soln_c[i][2][0], r_soln_c[i][2][1], r_soln_c[i][2][2]};
    solutions.push_back(sol);
  }

  return !solutions.empty();
}

} // namespace loop_closure
} // namespace f3dock
