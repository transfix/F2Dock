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
