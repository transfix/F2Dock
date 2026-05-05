#ifndef PROTEIN_FLEXIBILITY_LOOP_CLOSURE_TRIPEP_CLOSURE_H
#define PROTEIN_FLEXIBILITY_LOOP_CLOSURE_TRIPEP_CLOSURE_H

//!----------------------------------------------------------------------
//! Copyright (C) 2003 
//!      Chaok Seok, Evangelos Coutsias, Matthew Jacobson, and Ken Dill
//!      UCSF and Univeristy of New Mexico
//! Witten by Chaok Seok 2003.  
//! Ported to C++ by Nathan Clement 2017.
//!----------------------------------------------------------------------------
//!----------------------------------------------------------------------------
//!*************** Tripeptide Loop Closure Algorithm *****************
//! files to be compiled with:
//!    tripeptide_closure.f90
//!    sturm.c
//! 
//!*******************************************************************
//! subroutine  initialize_loop_closure(b_len, b_ang, t_ang)
//!*******************************************************************
//! connectivity of atoms:
//!   N1-A1-C1-N2-A2-C2-N3-A3-C3
//!
//! input:
//!
//  * b_len(1:6): bond lengths (A1-C1, C1-N2, ..., N3-A3)
//!  * b_ang(1:7): bond angles (N1-A1-C1, A1-C1-N2, ..., N3-A3-C3)
//!  * t_ang(1:2): torsion angles (A1-C1-N2-A2, A2-C2-N3-A3)
//!*******************************************************************
//!
//!*******************************************************************
//! subroutine solv_3pep_poly(r_n1, r_a1, r_a3, r_c3, &
//!     r_soln_n, r_soln_a, r_soln_c, n_soln)
//!*******************************************************************
//! input: 
//!  * r_n1(3), r_a1(3), r_a3(3), r_c3(3): 
//!       Cartesian coordinates of N and CA atoms of the first residue and
//!        CA and C atoms of the last (third) residue.
//! output:
//!  * n_soln: number of alternative loop closure solutions.
//!  * r_soln_n(3,3,8), r_soln_a(3,3,8), r_soln_c(3,3,8): 
//!       Cartesian coordinates of loop closure solutions. 
//!       first dimension: x, y, z component
//!       second dim: residue number
//!       third dim: solution number
//!*******************************************************************
//!----------------------------------------------------------------------------
//MODULE tripep_closure
//!----------------------------------------------------------------------------
//  integer, parameter :: dp = kind(1.0d0)
static constexpr double kPi = 3.141592653589793238462643383279502884197;
static constexpr double kTwoPi = 2.0 * kPi;
static constexpr double kDeg2Rad = kPi / 180.0;
static constexpr double kRad2Deg = 180.0 / kPi;
#define max_soln 16

class tripep_closure {
 public:
  //!-----------------------------------------------------------------------
  //subroutine initialize_loop_closure(b_len, b_ang, t_ang)
  void initialize_loop_closure(double b_len[6], double b_ang[7], double t_ang[2]);
  //!-----------------------------------------------------------------------
  //subroutine solv_3pep_poly(r_n1, r_a1, r_a3, r_c3, &
  //     r_soln_n, r_soln_a, r_soln_c, n_soln)
  void solve_3pep_poly(double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3], 
                       double r_soln_n[max_soln][3][3], 
                       double r_soln_a[max_soln][3][3], 
                       double r_soln_c[max_soln][3][3], int *n_soln);
  void set_print_level(int l) { print_level = l;};
  //!-----------------------------------------------------------------------

 private:
  //  integer, parameter :: print_level = 0
  int print_level = 0;
  //  ! parameters for tripeptide loop (including bond lengths & angles)
  //  real(dp) :: len0(6), b_ang0(7), t_ang0(2)
  double len0[6], b_ang0[7], t_ang0[2];
  //  real(dp) :: aa13_min_sqr, aa13_max_sqr
  double aa13_min_sqr, aa13_max_sqr;
  //  real(dp) :: delta(0:3), xi(3), eta(3), alpha(3), theta(3)
  double delta[4], xi[3], eta[3], alpha[3], theta[3];
  //  real(dp) :: cos_alpha(3), sin_alpha(3), cos_theta(3), sin_theta(3)
  double cos_alpha[3], sin_alpha[3], cos_theta[3], sin_theta[3];
  //  real(dp) :: cos_delta(0:3), sin_delta(0:3)
  double cos_delta[4], sin_delta[4];
  //  real(dp) :: cos_xi(3), cos_eta(3), sin_xi(3), sin_eta(3)
  double cos_xi[3], cos_eta[3], sin_xi[3], sin_eta[3];
  //  real(dp) :: r_a1a3(3), r_a1n1(3), r_a3c3(3)
  double r_a1a3[3], r_a1n1[3], r_a3c3[3];
  //  real(dp) :: b_a1a3(3), b_a1n1(3), b_a3c3(3)
  double b_a1a3[3], b_a1n1[3], b_a3c3[3];
  //  real(dp) :: len_na(3), len_ac(3), len_aa(3)
  double len_na[3], len_ac[3], len_aa[3];
  //  ! used for polynomial coefficients
  //  real(dp) :: C0(0:2,3), C1(0:2,3), C2(0:2,3)
  double C0[3][3], C1[3][3], C2[3][3];
  //  real(dp) :: Q(0:16,0:4), R(0:16,0:2)
  double Q[5][17], R[3][17];

  void get_input_angles(int *n_soln, 
                        double r_n1[3], double r_a1[3], 
                        double r_a3[3], double r_c3[3]);
  void test_two_cone_existence_soln(double tt, double kx, double et, double ap, 
                                    int *n_soln, char cone_type[2]);
  void get_poly_coeff(double poly_coeff[max_soln+1]);
  void poly_mul_sub2(double u1[5][5], double u2[5][5], 
                     double u3[5][5], double u4[5][5], 
                     int p1[2], int p2[2], 
                     int p3[2], int p4[2], 
                     double u5[5][5], int p5[2]);
  void poly_mul2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], 
                 double u3[5][5], int p3[2]);
  void poly_sub2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], 
                 double u3[5][5], int p3[2]);
  void poly_mul_sub1(double u1[17], double u2[17], double u3[17], double u4[17], 
                     int p1, int p2, int p3, int p4, double u5[17], int *p5);
  void poly_mul1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3);
  void poly_sub1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3);
  void coord_from_poly_roots(int *n_soln, double roots[max_soln], 
                             double r_n1[3], double r_a1[3], double r_a3[3], 
                             double r_c3[3], double r_soln_n[max_soln][3][3], 
                             double r_soln_a[max_soln][3][3], 
                             double r_soln_c[max_soln][3][3]);
  double calc_t2(double t0);
  double calc_t1(double t0, double t2);
  void calc_dih_ang(double r1[3], double r2[3], double r3[3], double *angle);
  void calc_bnd_ang(double r1[3], double r2[3], double *angle);
  void cross(double p[3], double q[3], double s[3]);
  void quaternion(double axis[3], double quarter_ang, double p[4]);
  void rotation_matrix(double q[4], double U[3][3]);
  //!----------------------------------------------------------------------------
  //END MODULE tripep_closure
  //!----------------------------------------------------------------------------
};
#endif // PROTEIN_FLEXIBILITY_LOOP_CLOSURE_TRIPEP_CLOSURE_H
