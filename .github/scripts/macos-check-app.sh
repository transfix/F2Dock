#!/usr/bin/env bash
# macos-check-app.sh APP [--run]
#
# Checks that a macOS .app made by macos-bundle-app.sh is self-contained.
#
# Always, statically: every library any executable in Contents/MacOS or any
# file in Contents/Frameworks loads is either part of the OS (/usr/lib,
# /System) or an @rpath/ name present in Contents/Frameworks, and every
# executable has the @executable_path/../Frameworks rpath.
#
# With --run (macOS only), also by launching each executable with nothing
# but the OS to fall back on: no DYLD_LIBRARY_PATH, DYLD_FALLBACK_LIBRARY_PATH
# limited to /usr/lib, and DYLD_PRINT_LIBRARIES on. Each must reach main()
# and print something, without a dyld error or a signal, and every image dyld
# loads must come from inside APP or from the OS. dyld binds every library
# an executable links before main() runs, so reaching main() proves the
# whole link-time closure resolves. Hide the build's dependency prefix
# first, or a library that is missing from the bundle may load from it.
#
# OTOOL overrides the tool, as in macos-bundle-app.sh.
#
# Kept bash-3.2-clean for macOS's /bin/bash.
set -euo pipefail

if [ $# -lt 1 ] || [ $# -gt 2 ] || { [ $# -eq 2 ] && [ "$2" != --run ]; }; then
  echo "usage: $0 APP [--run]" >&2
  exit 2
fi

OTOOL=${OTOOL:-otool}

abspath() { perl -MCwd -e 'print Cwd::abs_path(shift)' "$1"; }

app=$(abspath "$1")
exe_dir="$app/Contents/MacOS"
fw="$app/Contents/Frameworks"
bad=0

is_macho() { file -b "$1" | grep -q 'Mach-O'; }

refs() {
  local id
  id=$("$OTOOL" -D "$1" | sed -n '2p')
  "$OTOOL" -L "$1" | sed '1d' | awk '{print $1}' | while read -r r; do
    if [ "$r" != "$id" ]; then echo "$r"; fi
  done
}

rpaths() {
  "$OTOOL" -l "$1" | awk '$1 == "cmd" && $2 == "LC_RPATH" { want = 1; next }
    want && $1 == "path" { print $2; want = 0 }'
}

# Frameworks/ is walked recursively: a framework keeps its binary under
# Versions/<v>/.
files=$(
  find "$exe_dir" "$fw" -type f | sort | while read -r f; do
    if is_macho "$f"; then echo "$f"; fi
  done
)

for f in $files; do
  for r in $(refs "$f"); do
    case "$r" in
      /usr/lib/* | /System/*) ;;
      @rpath/*)
        if [ ! -e "$fw/${r#@rpath/}" ]; then
          echo "FAIL: ${f#"$app"/} loads $r, which is not in Contents/Frameworks"
          bad=1
        fi
        ;;
      *)
        echo "FAIL: ${f#"$app"/} loads $r, which is neither the OS's nor @rpath/"
        bad=1
        ;;
    esac
  done
done

exes=$(
  for f in "$exe_dir"/*; do
    if [ -f "$f" ] && [ ! -L "$f" ] && is_macho "$f"; then echo "$f"; fi
  done
)
for f in $exes; do
  if ! rpaths "$f" | grep -qxF '@executable_path/../Frameworks'; then
    echo "FAIL: ${f#"$app"/} has no @executable_path/../Frameworks rpath"
    bad=1
  fi
done

if [ "$bad" -ne 0 ]; then exit 1; fi
echo "static: $(echo "$files" | wc -l | tr -d ' ') Mach-O files, every reference resolves inside the bundle or the OS"

if [ $# -lt 2 ]; then exit 0; fi

if [ "$(uname -s)" != Darwin ]; then
  echo "--run needs macOS" >&2
  exit 2
fi

out=$(mktemp)
err=$(mktemp)
trap 'rm -f "$out" "$err"' EXIT

for f in $exes; do
  name=$(basename "$f")
  # F2Dock has a --version; the others print their usage when given nothing.
  # None of them starts work, or a server, in that case.
  args=""
  if [ "$name" = F2Dock ]; then args=--version; fi
  # perl is the timeout and must come first: /usr/bin/perl and /usr/bin/env
  # are SIP-protected, and macOS strips DYLD_* from the environment of a
  # protected binary it launches, so env has to be the last thing to set
  # them before the binary under test. The alarm survives both execs.
  status=0
  perl -e 'alarm 60; exec @ARGV or die "exec: $!"' \
    env -u DYLD_LIBRARY_PATH -u DYLD_FRAMEWORK_PATH -u DYLD_INSERT_LIBRARIES \
    DYLD_FALLBACK_LIBRARY_PATH=/usr/lib DYLD_PRINT_LIBRARIES=1 \
    "$f" $args >"$out" 2>"$err" || status=$?

  loaded=$(sed -n 's/^dyld\[[0-9]*\]: <[^>]*> //p' "$err")
  outside=$(echo "$loaded" | awk -v a="$app/" 'NF && index($0, a) != 1 && $0 !~ /^\/usr\/lib\// && $0 !~ /^\/System\//')
  said=$({ cat "$out"; grep -v '^dyld\[' "$err" || true; } | sed -n '1p')
  ok=yes
  if [ "$status" -ge 128 ]; then
    echo "FAIL: $name died with signal $((status - 128))"
    ok=no
  fi
  if grep -q 'Library not loaded\|Symbol not found\|dyld\[[0-9]*\]: [A-Z]' "$err"; then
    echo "FAIL: $name: dyld error:"
    grep -v '^dyld\[[0-9]*\]: <' "$err" | head -n 5
    ok=no
  fi
  if [ -n "$outside" ]; then
    echo "FAIL: $name loaded images from outside the bundle and the OS:"
    echo "$outside"
    ok=no
  fi
  if [ -z "$said" ]; then
    echo "FAIL: $name printed nothing"
    ok=no
  fi
  if [ "$ok" = yes ]; then
    n=$(echo "$loaded" | awk -v a="$app/" 'index($0, a) == 1' | wc -l | tr -d ' ')
    echo "run: $name $args -> exit $status, $n images from the bundle: $said"
  else
    bad=1
  fi
done

exit "$bad"
