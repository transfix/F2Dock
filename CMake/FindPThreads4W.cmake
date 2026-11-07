#
# FindPThreads4W -- pthreads4w (POSIX threads for Windows) as the
# PThreads4W::PThreads4W target the F2Dock targets link on Windows.
#
# vcpkg's port shipped a PThreads4W package config, which is where that target
# name comes from. The cvcpkg `pthreads4w` bundle ships none -- only
# include/pthread.h, lib/pthreadVC3.lib and bin/pthreadVC3.dll -- so a config
# is used when one is on the search path, and otherwise the library is found
# directly and the same target is defined over it.
#
# The include directory is located through _ptw32.h, which only pthreads4w
# ships, so it cannot resolve to the MinGW winpthreads pthread.h that other
# bundles in the same prefix may carry.
#
# Result: PThreads4W_FOUND and the PThreads4W::PThreads4W target.
#

find_package(PThreads4W CONFIG QUIET)

if(NOT TARGET PThreads4W::PThreads4W)
  find_path(PThreads4W_INCLUDE_DIR NAMES _ptw32.h)
  find_library(PThreads4W_LIBRARY NAMES pthreadVC3 pthreadVC3d)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(PThreads4W
    REQUIRED_VARS PThreads4W_LIBRARY PThreads4W_INCLUDE_DIR)

  if(PThreads4W_FOUND)
    add_library(PThreads4W::PThreads4W UNKNOWN IMPORTED)
    set_target_properties(PThreads4W::PThreads4W PROPERTIES
      IMPORTED_LOCATION "${PThreads4W_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${PThreads4W_INCLUDE_DIR}")
  endif()

  mark_as_advanced(PThreads4W_INCLUDE_DIR PThreads4W_LIBRARY)
else()
  set(PThreads4W_FOUND TRUE)
endif()
