/*
  Copyright 2011 The University of Texas at Austin

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
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
  USA
*/

// RAWIV.cpp now delegates volume serialization to libcvc's volume I/O stack
// (see cvc/volume.h, cvc/volume_file_info.h). libcvc handles RAWIV header
// construction, endian swapping, and per-voxel data-type packing for UChar /
// Float / Double, matching this module's prior behavior. The public API
// (writeGrid / readShapeCompGrid / readElecGrid / readRAWIVHeader) is
// preserved unchanged so existing F2Dock callers do not need to be updated.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cstddef>
#include <exception>
#include <string>
#include <vector>

#include <cvc/app.h>
#include <cvc/bounding_box.h>
#include <cvc/dimension.h>
#include <cvc/types.h>
#include <cvc/volume.h>
#include <cvc/volume_file_info.h>

#include "RAWIV.h"

namespace {

#ifdef FFTW_SINGLE_PRECISION
constexpr CVC_NAMESPACE::data_type kFftwVoxelType = CVC_NAMESPACE::Float;
#else
constexpr CVC_NAMESPACE::data_type kFftwVoxelType = CVC_NAMESPACE::Double;
#endif

// Build a libcvc bounding_box from the (xCenter, yCenter, zCenter, scale)
// convention used throughout F2Dock. F2Dock stores `scale` as the inverse of
// the cube edge length, so the half-extent in each dimension is 0.5/scale.
CVC_NAMESPACE::bounding_box centerScaleToBox(double xCenter, double yCenter,
                                             double zCenter, double scale) {
  const double half = 0.5 / scale;
  return CVC_NAMESPACE::bounding_box(xCenter - half, yCenter - half,
                                     zCenter - half, xCenter + half,
                                     yCenter + half, zCenter + half);
}

// Helper: write an n^3 cube of FFTW samples to a RAWIV file via libcvc.
int writeRAWIVCube(const FFTW_DATA_TYPE *vol, int n, double xCenter,
                   double yCenter, double zCenter, double scale,
                   const char *fileName) {
  if (vol == NULL || n < 1 || scale <= 0 || fileName == NULL) {
    return 0;
  }

  try {
    CVC_NAMESPACE::app ctx;
    CVC_NAMESPACE::volume cvol(
        ctx, reinterpret_cast<const unsigned char *>(vol),
        CVC_NAMESPACE::dimension(n, n, n), kFftwVoxelType,
        centerScaleToBox(xCenter, yCenter, zCenter, scale));
    cvol.write(std::string(fileName));
    return 1;
  } catch (const std::exception &e) {
    fprintf(stderr, "\nError: Failed to write RAWIV file %s (%s)!\n\n",
            fileName, e.what());
    return 0;
  }
}

// Helper: load a RAWIV file via libcvc and copy the voxel data into a
// freshly-allocated FFTW_DATA_TYPE buffer (caller takes ownership via
// free()). Returns the cube dimensions and (center, scale) recovered from
// the file's bounding box.
int readRAWIVCube(FFTW_DATA_TYPE **outVol, int *xDim, int *yDim, int *zDim,
                  double *xCenter, double *yCenter, double *zCenter,
                  double *scale, const char *fileName) {
  if (outVol == NULL || fileName == NULL) {
    return 0;
  }

  try {
    CVC_NAMESPACE::app ctx;
    CVC_NAMESPACE::volume cvol(ctx, std::string(fileName));

    *xDim = static_cast<int>(cvol.XDim());
    *yDim = static_cast<int>(cvol.YDim());
    *zDim = static_cast<int>(cvol.ZDim());

    *xCenter = (cvol.XMin() + cvol.XMax()) / 2.0;
    *yCenter = (cvol.YMin() + cvol.YMax()) / 2.0;
    *zCenter = (cvol.ZMin() + cvol.ZMax()) / 2.0;
    *scale = 1.0 / (cvol.XMax() - cvol.XMin());

    const std::size_t n = static_cast<std::size_t>(*xDim) *
                          static_cast<std::size_t>(*yDim) *
                          static_cast<std::size_t>(*zDim);

    *outVol = static_cast<FFTW_DATA_TYPE *>(malloc(n * sizeof(FFTW_DATA_TYPE)));
    if (*outVol == NULL) {
      fprintf(
          stderr,
          "\nError: Failed to allocate buffer space for reading file %s!\n\n",
          fileName);
      return 0;
    }

    // libcvc converts on read to whatever voxelType() the file declared.
    // Cast through the per-element accessor so we get a clean conversion to
    // FFTW_DATA_TYPE regardless of whether the on-disk type was Float,
    // Double, or UChar.
    const CVC_NAMESPACE::data_type vt = cvol.voxelType();
    const unsigned char *raw = *cvol;

    if (vt == kFftwVoxelType) {
      memcpy(*outVol, raw, n * sizeof(FFTW_DATA_TYPE));
    } else if (vt == CVC_NAMESPACE::Float) {
      const float *src = reinterpret_cast<const float *>(raw);
      for (std::size_t i = 0; i < n; ++i)
        (*outVol)[i] = static_cast<FFTW_DATA_TYPE>(src[i]);
    } else if (vt == CVC_NAMESPACE::Double) {
      const double *src = reinterpret_cast<const double *>(raw);
      for (std::size_t i = 0; i < n; ++i)
        (*outVol)[i] = static_cast<FFTW_DATA_TYPE>(src[i]);
    } else if (vt == CVC_NAMESPACE::UChar) {
      const unsigned char *src = raw;
      for (std::size_t i = 0; i < n; ++i)
        (*outVol)[i] = static_cast<FFTW_DATA_TYPE>(src[i]);
    } else {
      fprintf(
          stderr,
          "\nError: Unsupported voxel type in file %s for FFTW grid load!\n\n",
          fileName);
      free(*outVol);
      *outVol = NULL;
      return 0;
    }

    printf("\n\nREAD %s ( xDim = %d, yDim = %d, zDim = %d, xCenter = %lf, "
           "yCenter = %lf, zCenter = %lf, scale = %lf, voxelType = %s )!\n\n",
           fileName, *xDim, *yDim, *zDim, *xCenter, *yCenter, *zCenter, *scale,
           cvol.voxelTypeStr());
    return 1;
  } catch (const std::exception &e) {
    fprintf(stderr, "\nError: Failed to read RAWIV file %s (%s)!\n\n", fileName,
            e.what());
    return 0;
  }
}

} // namespace

int writeGrid(FFTW_complex *scGrid, FFTW_DATA_TYPE *elecGrid, int n,
              double xCenter, double yCenter, double zCenter, double scale,
              char *fileName, char *fileNameSCRe, char *fileNameSCIm,
              char *fileNameElecRe) {
  if ((scGrid == NULL) || (elecGrid == NULL) || (n < 1) || (scale <= 0) ||
      (fileName == NULL))
    return 0;

  const int nCube = n * n * n;

  // Split scGrid into separate real and imaginary RAWIV files. Use plain
  // std::vector buffers so libcvc copies them into its own storage.
  std::vector<FFTW_DATA_TYPE> reBuf(nCube);
  std::vector<FFTW_DATA_TYPE> imBuf(nCube);
  for (int i = 0; i < nCube; i++) {
    reBuf[i] = scGrid[i][0];
    imBuf[i] = scGrid[i][1];
  }

  // Build default per-channel filenames (<fileName>-SC-Re.rawiv, etc.) when
  // the caller passes NULL for any of the three. Mirrors the prior behavior
  // exactly so callers (Docking.cpp) keep working.
  const std::size_t l = strlen(fileName);
  const std::size_t bufLen = l + 15;
  std::vector<char> bufSCRe(bufLen), bufSCIm(bufLen), bufElecRe(bufLen);
  char *fileNameSCRe2 = bufSCRe.data();
  char *fileNameSCIm2 = bufSCIm.data();
  char *fileNameElecRe2 = bufElecRe.data();

  strcpy(fileNameSCRe2, fileName);
  strcpy(fileNameSCIm2, fileName);
  strcpy(fileNameElecRe2, fileName);

  if ((l > 3) && (fileName[l - 4] == '.'))
    fileNameSCRe2[l - 4] = fileNameSCIm2[l - 4] = fileNameElecRe2[l - 4] = 0;

  strcat(fileNameSCRe2, "-SC-Re.rawiv");
  strcat(fileNameSCIm2, "-SC-Im.rawiv");
  strcat(fileNameElecRe2, "-Elec-Re.rawiv");

  if (fileNameSCRe == NULL)
    fileNameSCRe = fileNameSCRe2;
  if (fileNameSCIm == NULL)
    fileNameSCIm = fileNameSCIm2;
  if (fileNameElecRe == NULL)
    fileNameElecRe = fileNameElecRe2;

  if (!writeRAWIVCube(reBuf.data(), n, xCenter, yCenter, zCenter, scale,
                      fileNameSCRe) ||
      !writeRAWIVCube(imBuf.data(), n, xCenter, yCenter, zCenter, scale,
                      fileNameSCIm) ||
      !writeRAWIVCube(elecGrid, n, xCenter, yCenter, zCenter, scale,
                      fileNameElecRe)) {
    return 0;
  }

  return 1;
}

int readRAWIVHeader(int *xDim, int *yDim, int *zDim, double *xCenter,
                    double *yCenter, double *zCenter, double *scale,
                    char *fileName) {
  if (fileName == NULL)
    return 0;

  try {
    CVC_NAMESPACE::app ctx;
    CVC_NAMESPACE::volume_file_info info(ctx, std::string(fileName));

    *xDim = static_cast<int>(info.XDim());
    *yDim = static_cast<int>(info.YDim());
    *zDim = static_cast<int>(info.ZDim());

    *xCenter = (info.boundingBox().minx + info.boundingBox().maxx) / 2.0;
    *yCenter = (info.boundingBox().miny + info.boundingBox().maxy) / 2.0;
    *zCenter = (info.boundingBox().minz + info.boundingBox().maxz) / 2.0;
    *scale = 1.0 / (info.boundingBox().maxx - info.boundingBox().minx);

    return 1;
  } catch (const std::exception &e) {
    fprintf(stderr, "\nError: Failed to read header from file %s (%s)!\n\n",
            fileName, e.what());
    return 0;
  }
}

int readShapeCompGrid(FFTW_complex **scGrid, int *n, double *xCenter,
                      double *yCenter, double *zCenter, double *scale,
                      char *fileNameRe, char *fileNameIm) {
  if ((fileNameRe == NULL) || (fileNameIm == NULL))
    return 0;

  FFTW_DATA_TYPE *vol = NULL;
  int xDim, yDim, zDim;

  if (!readRAWIVCube(&vol, &xDim, &yDim, &zDim, xCenter, yCenter, zCenter,
                     scale, fileNameRe))
    return 0;

  if ((xDim != yDim) || (yDim != zDim) || (zDim != xDim)) {
    free(vol);
    return 0;
  }

  *n = xDim;
  const std::size_t n3 = static_cast<std::size_t>(*n) *
                         static_cast<std::size_t>(*n) *
                         static_cast<std::size_t>(*n);
  *scGrid = static_cast<FFTW_complex *>(FFTW_malloc(n3 * sizeof(FFTW_complex)));

  if (*scGrid == NULL) {
    free(vol);
    return 0;
  }

  for (std::size_t c = 0; c < n3; c++)
    (*scGrid)[c][0] = vol[c];

  free(vol);
  vol = NULL;

  double xc, yc, zc, sc;

  if (!readRAWIVCube(&vol, &xDim, &yDim, &zDim, &xc, &yc, &zc, &sc, fileNameIm))
    return 0;

  if ((xDim != (*n)) || (yDim != (*n)) || (zDim != (*n)) ||
      (xc != (*xCenter)) || (yc != (*yCenter)) || (zc != (*zCenter)) ||
      (sc != (*scale))) {
    printf("\nError: Parameter mismatch between files %s and %s!\n", fileNameRe,
           fileNameIm);
    FFTW_free(*scGrid);
    *scGrid = NULL;
    free(vol);
    return 0;
  }

  for (std::size_t c = 0; c < n3; c++)
    (*scGrid)[c][1] = vol[c];

  free(vol);

  return 1;
}

int readShapeCompGrid(FFTW_complex *scGrid, double *xCenter, double *yCenter,
                      double *zCenter, char *fileNameRe, char *fileNameIm) {
  if ((scGrid == NULL) || (fileNameRe == NULL) || (fileNameIm == NULL))
    return 0;

  FFTW_DATA_TYPE *vol = NULL;
  int xDim, yDim, zDim;
  double scale;

  if (!readRAWIVCube(&vol, &xDim, &yDim, &zDim, xCenter, yCenter, zCenter,
                     &scale, fileNameRe))
    return 0;

  if ((xDim != yDim) || (yDim != zDim) || (zDim != xDim)) {
    free(vol);
    return 0;
  }

  const int n = xDim;
  const std::size_t n3 = static_cast<std::size_t>(n) *
                         static_cast<std::size_t>(n) *
                         static_cast<std::size_t>(n);

  for (std::size_t c = 0; c < n3; c++)
    scGrid[c][0] = vol[c];

  free(vol);
  vol = NULL;

  double xc, yc, zc, sc;

  if (!readRAWIVCube(&vol, &xDim, &yDim, &zDim, &xc, &yc, &zc, &sc, fileNameIm))
    return 0;

  if ((xDim != n) || (yDim != n) || (zDim != n) || (xc != (*xCenter)) ||
      (yc != (*yCenter)) || (zc != (*zCenter)) || (sc != scale)) {
    printf("\nError: Parameter mismatch between files %s and %s!\n", fileNameRe,
           fileNameIm);
    free(vol);
    return 0;
  }

  for (std::size_t c = 0; c < n3; c++)
    scGrid[c][1] = vol[c];

  free(vol);

  return 1;
}

int readElecGrid(FFTW_DATA_TYPE **elecGrid, int *n, double *xCenter,
                 double *yCenter, double *zCenter, double *scale,
                 char *fileNameRe) {
  if (fileNameRe == NULL)
    return 0;

  FFTW_DATA_TYPE *vol = NULL;
  int xDim, yDim, zDim;

  if (!readRAWIVCube(&vol, &xDim, &yDim, &zDim, xCenter, yCenter, zCenter,
                     scale, fileNameRe))
    return 0;

  if ((xDim != yDim) || (yDim != zDim) || (zDim != xDim)) {
    free(vol);
    return 0;
  }

  *n = xDim;
  const std::size_t n3 = static_cast<std::size_t>(*n) *
                         static_cast<std::size_t>(*n) *
                         static_cast<std::size_t>(*n);
  // Preserve the original sizing quirk: the prior implementation allocated
  // n^3 * sizeof(FFTW_complex) here even though only the real part is used.
  // Match that allocation size to keep callers' assumptions intact.
  *elecGrid =
      static_cast<FFTW_DATA_TYPE *>(FFTW_malloc(n3 * sizeof(FFTW_complex)));

  if (*elecGrid == NULL) {
    free(vol);
    return 0;
  }

  for (std::size_t c = 0; c < n3; c++)
    (*elecGrid)[c] = vol[c];

  free(vol);

  return 1;
}

int readElecGrid(FFTW_DATA_TYPE *elecGrid, double *xCenter, double *yCenter,
                 double *zCenter, char *fileNameRe) {
  if ((elecGrid == NULL) || (fileNameRe == NULL))
    return 0;

  FFTW_DATA_TYPE *vol = NULL;
  int xDim, yDim, zDim;
  double scale;

  if (!readRAWIVCube(&vol, &xDim, &yDim, &zDim, xCenter, yCenter, zCenter,
                     &scale, fileNameRe))
    return 0;

  if ((xDim != yDim) || (yDim != zDim) || (zDim != xDim)) {
    free(vol);
    return 0;
  }

  const int n = xDim;
  const std::size_t n3 = static_cast<std::size_t>(n) *
                         static_cast<std::size_t>(n) *
                         static_cast<std::size_t>(n);

  for (std::size_t c = 0; c < n3; c++)
    elecGrid[c] = vol[c];

  free(vol);

  return 1;
}
