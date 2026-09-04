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

/*
  The FFTW_* names below are the engine's entire FFT vocabulary -- Docking.cpp
  and fastfft.cpp never call a transform library directly. That makes this the
  one place a provider is chosen, and why swapping FFTW out costs no churn in
  the 8000-line scoring loop.

  F2DOCK_FFT_PROVIDER_FFTW selects FFTW itself and maps these straight onto it,
  as before. It is OPT-IN: FFTW is GPL-2.0-or-later, so that build produces a
  GPL binary rather than an LGPL-2.1 one. Every other provider goes through
  cvc::fft (see fft_provider.h), which deliberately keeps FFTW's semantics --
  unnormalised both ways, half-complex r2c layout, plans that capture their
  buffers -- so the choice of provider is not a change in behaviour.
*/

#ifdef F2DOCK_FFT_PROVIDER_FFTW

#include <fftw3.h>

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

#else /* provider-neutral path: PocketFFT, MKL, cuFFT, ... */

#include "fft-utils/fft_provider.h"

#define FFTW_DATA_TYPE cvc::fft::real
#define FFTW_complex cvc::fft::complex
#define FFTW_plan cvc::fft::plan
#define FFTW_execute cvc::fft::execute
#define FFTW_free cvc::fft::release
#define FFTW_malloc cvc::fft::alloc
#define FFTW_plan_dft_1d cvc::fft::dft_1d
#define FFTW_plan_dft_2d cvc::fft::dft_2d
#define FFTW_plan_dft_3d cvc::fft::dft_3d
#define FFTW_plan_dft_r2c_3d cvc::fft::dft_r2c_3d
#define FFTW_plan_dft_c2r_3d cvc::fft::dft_c2r_3d
#define FFTW_destroy_plan cvc::fft::destroy

/*
  Deliberately NOT defined here: FFTW_plan_many_dft, FFTW_execute_dft and the
  wisdom pair. The first two were used only by the GPL sparsefft3 module and
  the wisdom pair only by rank-fftw's grid-timing utility, both of which are
  FFTW-only. Leaving them undefined means any future use outside an FFTW build
  is a compile error naming the symbol, rather than a silent fallback.
*/

/* FFTW's own numeric values, so a plan flag or sign written for one provider
   means the same thing under another. */
#ifndef FFTW_FORWARD
#define FFTW_FORWARD (-1)
#endif
#ifndef FFTW_BACKWARD
#define FFTW_BACKWARD (1)
#endif
#ifndef FFTW_MEASURE
#define FFTW_MEASURE (0U)
#endif
#ifndef FFTW_ESTIMATE
#define FFTW_ESTIMATE (1U << 6)
#endif
#ifndef FFTW_PATIENT
#define FFTW_PATIENT (1U << 5)
#endif
/* rank-fftw passes this when timing a size; providers with no planning stage
   ignore it, and nothing reads the input array afterwards either way. */
#ifndef FFTW_DESTROY_INPUT
#define FFTW_DESTROY_INPUT (1U << 0)
#endif

#endif /* F2DOCK_FFT_PROVIDER_FFTW */

#ifndef OPT_FFTW_SEARCH_TYPE
#define OPT_FFTW_SEARCH_TYPE FFTW_MEASURE
#endif
#endif
