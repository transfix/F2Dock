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

#ifndef F3DOCK_UTIL_STRING_UTIL_H
#define F3DOCK_UTIL_STRING_UTIL_H

#include <algorithm>
#include <cctype>
#include <string_view>

namespace f3dock {
namespace util {

// Case-insensitive ASCII string equality. Portable replacement for
// POSIX strcasecmp / MSVC _stricmp; relies only on the C++20 standard
// library.
inline bool iequals(std::string_view a, std::string_view b) {
  return a.size() == b.size() &&
         std::equal(a.begin(), a.end(), b.begin(),
                    [](unsigned char ac, unsigned char bc) {
                      return std::tolower(ac) == std::tolower(bc);
                    });
}

} // namespace util
} // namespace f3dock

#endif // F3DOCK_UTIL_STRING_UTIL_H
