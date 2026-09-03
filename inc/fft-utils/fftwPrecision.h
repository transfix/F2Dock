/*
  Copyright 2011 The University of Texas at Austin

        Authors: Rezaul Alam Chowdhury <shaikat@cs.utexas.edu>
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

#ifndef FFTW_PRECISION_DEFINED
#define FFTW_PRECISION_DEFINED

#ifdef FFTW_SINGLE_PRECISION
#define FFTW_DATA_TYPE float
#define FFTW_complex fftwf_complex
#define FFTW_plan fftwf_plan
#define FFTW_execute fftwf_execute
#define FFTW_execute_dft fftwf_execute_dft
#define FFTW_free fftwf_free
#define FFTW_malloc fftwf_malloc
#define FFTW_plan_many_dft fftwf_plan_many_dft
#define FFTW_plan_dft_1d fftwf_plan_dft_1d
#define FFTW_plan_dft_2d fftwf_plan_dft_2d
#define FFTW_plan_dft_3d fftwf_plan_dft_3d
#define FFTW_plan_dft_r2c_3d fftwf_plan_dft_r2c_3d
#define FFTW_plan_dft_c2r_3d fftwf_plan_dft_c2r_3d
#define FFTW_destroy_plan fftwf_destroy_plan
#define FFTW_export_wisdom_to_file fftwf_export_wisdom_to_file
#define FFTW_import_wisdom_from_file fftwf_import_wisdom_from_file
#else
#define FFTW_DATA_TYPE double
#define FFTW_complex fftw_complex
#define FFTW_plan fftw_plan
#define FFTW_execute fftw_execute
#define FFTW_execute_dft fftw_execute_dft
#define FFTW_free fftw_free
#define FFTW_malloc fftw_malloc
#define FFTW_plan_many_dft fftw_plan_many_dft
#define FFTW_plan_dft_1d fftw_plan_dft_1d
#define FFTW_plan_dft_2d fftw_plan_dft_2d
#define FFTW_plan_dft_3d fftw_plan_dft_3d
#define FFTW_plan_dft_r2c_3d fftw_plan_dft_r2c_3d
#define FFTW_plan_dft_c2r_3d fftw_plan_dft_c2r_3d
#define FFTW_destroy_plan fftw_destroy_plan
#define FFTW_export_wisdom_to_file fftw_export_wisdom_to_file
#define FFTW_import_wisdom_from_file fftw_import_wisdom_from_file
#endif

#ifndef OPT_FFTW_SEARCH_TYPE
#define OPT_FFTW_SEARCH_TYPE FFTW_MEASURE
#endif
#endif
