# cvcpkg/recipes/f2dock/build.ps1 — build F2Dock against an installed libcvc SDK.
#
# Windows goes through _common/env-windows.ps1 rather than a hand-rolled
# `-G Ninja`: it forces CC/CXX=cl, imports the MSVC developer environment,
# strips MinGW/MSYS dirs from PATH, and sets CMAKE_MSVC_RUNTIME_LIBRARY so the
# CRT matches the link mode. That is an ABI requirement -- every cvcpkg windows
# dep bundle is MSVC-built, so a gcc-compiled object could not link them, and
# CMake would otherwise pick C:\mingw64\bin\c++.exe off the runner PATH.
$ErrorActionPreference = 'Stop'

if (-not $env:CVC_SOURCE_DIR) { throw 'CVC_SOURCE_DIR must be set' }
if (-not $env:CVC_BUILD_DIR) { throw 'CVC_BUILD_DIR must be set' }
if (-not $env:CVC_INSTALL_DIR) { throw 'CVC_INSTALL_DIR must be set' }

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
. "$scriptDir\..\_common\env-windows.ps1"

# CMAKE_PREFIX_PATH is deliberately NOT passed: env-windows.ps1 sets it to
# CVC_DEPS_PREFIX;CVC_BUILD_PREFIX, and a -D here would override that and drop
# the build prefix.
Invoke-CvcCMakeBuild @(
    '-DF2DOCK_BUILD_TESTS=OFF'
    '-DF2DOCK_ENABLE_GPL_COMPONENTS=OFF'
)
