cmake_minimum_required(VERSION 3.16)

ADD_SUBDIRECTORY (fft-utils)
ADD_SUBDIRECTORY (math)
ADD_SUBDIRECTORY (misc-ident)
ADD_SUBDIRECTORY (utils)
ADD_SUBDIRECTORY (vol)
ADD_SUBDIRECTORY (fast-clash)
ADD_SUBDIRECTORY (fast-GB)
ADD_SUBDIRECTORY (fast-hydro)
ADD_SUBDIRECTORY (fast-LJ)
ADD_SUBDIRECTORY (fast-PQ)
ADD_SUBDIRECTORY (fast-resCont)

OPTION(HAVE_LIBMOL "Whether or not the libMol package is part of this build" OFF)
IF(HAVE_LIBMOL)
ADD_DEFINITIONS (-DLIBMOL_FOUND)
ADD_SUBDIRECTORY (libmol)
ADD_SUBDIRECTORY (hbondFilter)
ENDIF(HAVE_LIBMOL)

ADD_SUBDIRECTORY (PG-range)
if(F2DOCK_ENABLE_F3DOCK_EXPERIMENTAL)
ADD_SUBDIRECTORY (f3dock-loop-closure)
# libicp (src/f3dock-icp/{icp,icpPointToPoint,icpPointToPlane,matrix}.cpp
# and inc/libicp/) is GPL-2.0-or-later, and kdtree.cpp is AFL-1.1. The whole
# target is gated: IcpAligner is a thin wrapper over libicp and has nothing
# to wrap without it.
if(F2DOCK_ENABLE_GPL_COMPONENTS)
ADD_SUBDIRECTORY (f3dock-icp)
endif()
ADD_SUBDIRECTORY (f3dock-flex)
ADD_SUBDIRECTORY (f3dock-domain)
ADD_SUBDIRECTORY (f3dock-nfft)
endif()
ADD_SUBDIRECTORY (f2dock)
ADD_SUBDIRECTORY (GB-rerank)

if(TARGET cvc::xmlrpc)
  ADD_SUBDIRECTORY (F2DockServer)
else()
  message(FATAL_ERROR "cvc::xmlrpc target not found; libcvc XML-RPC dependency was not configured")
endif()
ADD_SUBDIRECTORY (PostFiltering)
