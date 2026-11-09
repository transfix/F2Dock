#!/usr/bin/env bash
# macos-bundle-app.sh APP PREFIX
#
# Makes a macOS .app self-contained. Every dylib its executables load that is
# not part of the OS is copied into APP/Contents/Frameworks, transitively, and
# each executable gets an LC_RPATH of @executable_path/../Frameworks.
#
# PREFIX is the cvcpkg prefix the binaries were built against. A non-system
# library that resolves anywhere else is an error, not something to bundle:
# it would mean the build picked up a library from outside the hermetic
# prefix (Homebrew, a stray /usr/local).
#
# Libraries are copied under their real (versioned) names, and every other
# name they are referenced by becomes a symlink to that one copy, as in the
# prefix: libz.1.dylib -> libz.1.3.1.dylib. A library reached under two
# names is therefore still one file and one loaded image. A framework
# (@rpath/png.framework/Versions/1.6.44/png) is copied whole, symlinks and
# signature included, so its reference resolves unchanged.
#
# Every Mach-O file this changes is re-signed ad hoc. install_name_tool
# invalidates the signature, and arm64 macOS kills a binary whose signature
# does not verify.
#
# OTOOL, INSTALL_NAME_TOOL and CODESIGN override the tools, so the walk can
# be exercised off macOS (llvm-otool, llvm-install-name-tool, CODESIGN=true).
#
# Kept bash-3.2-clean for macOS's /bin/bash.
set -euo pipefail

if [ $# -ne 2 ]; then
  echo "usage: $0 APP PREFIX" >&2
  exit 2
fi

OTOOL=${OTOOL:-otool}
INSTALL_NAME_TOOL=${INSTALL_NAME_TOOL:-install_name_tool}
CODESIGN=${CODESIGN:-codesign}

# Resolved, symlink-free absolute path (macOS 12 has no realpath(1)).
abspath() { perl -MCwd -e 'print Cwd::abs_path(shift)' "$1"; }

app=$(abspath "$1")
prefix=$(abspath "$2")
exe_dir="$app/Contents/MacOS"
fw="$app/Contents/Frameworks"
mkdir -p "$fw"

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
# One line per Mach-O file still to scan: "<file><TAB><dir its @loader_path
# meant where it came from>".
: >"$work/queue"
# Files whose load commands changed, to re-sign.
: >"$work/modified"

is_macho() { file -b "$1" | grep -q 'Mach-O'; }

is_system() {
  case "$1" in
    /usr/lib/* | /System/*) return 0 ;;
    *) return 1 ;;
  esac
}

# The libraries FILE loads, one per line: otool -L's list minus its header
# line and, for a dylib, minus its own install name.
refs() {
  local id
  id=$("$OTOOL" -D "$1" | sed -n '2p')
  "$OTOOL" -L "$1" | sed '1d' | awk '{print $1}' | while read -r r; do
    if [ "$r" != "$id" ]; then echo "$r"; fi
  done
}

# FILE's LC_RPATH entries, one per line.
rpaths() {
  "$OTOOL" -l "$1" | awk '$1 == "cmd" && $2 == "LC_RPATH" { want = 1; next }
    want && $1 == "path" { print $2; want = 0 }'
}

# resolve REF FILE SRCDIR: the path in the prefix that REF, as loaded by
# FILE (which came from SRCDIR), refers to. Empty if nothing matches.
resolve() {
  local r=$1 f=$2 src=$3 rest p cand
  case "$r" in
    @rpath/*)
      rest=${r#@rpath/}
      for p in $(rpaths "$f") "$prefix/lib"; do
        case "$p" in
          @loader_path*) p="$src${p#@loader_path}" ;;
          @executable_path*) continue ;;
        esac
        cand="$p/$rest"
        if [ -e "$cand" ]; then
          echo "$cand"
          return 0
        fi
      done
      ;;
    @loader_path/*)
      cand="$src/${r#@loader_path/}"
      if [ -e "$cand" ]; then echo "$cand"; fi
      ;;
    /*)
      if [ -e "$r" ]; then echo "$r"; fi
      ;;
  esac
}

in_prefix() {
  case "$1" in
    "$prefix"/*) ;;
    *)
      echo "error: $2 resolves to $1, outside the cvcpkg prefix $prefix" >&2
      exit 1
      ;;
  esac
}

# bundle_framework PATH REST: PATH is where @rpath/REST resolved, and REST
# names a file inside X.framework. Copy the whole framework, once.
bundle_framework() {
  local path=$1 rest=$2 fwrel base
  fwrel=${rest%%.framework/*}.framework
  base=${path%"/$rest"}
  in_prefix "$(abspath "$base/$fwrel")" "@rpath/$rest"
  if [ ! -e "$fw/$fwrel" ]; then
    mkdir -p "$(dirname "$fw/$fwrel")"
    cp -RP "$base/$fwrel" "$fw/$fwrel"
    printf '%s\t%s\n' "$(abspath "$fw/$rest")" "$(dirname "$(abspath "$path")")" >>"$work/queue"
  fi
}

# bundle SRC NAME: copy the library SRC resolves to into Frameworks under its
# real name, once, and make NAME a symlink to it if NAME differs.
bundle() {
  local real name=$2 realname
  real=$(abspath "$1")
  in_prefix "$real" "$name"
  realname=$(basename "$real")
  if [ ! -e "$fw/$realname" ]; then
    cp "$real" "$fw/$realname"
    chmod u+w "$fw/$realname"
    printf '%s\t%s\n' "$fw/$realname" "$(dirname "$real")" >>"$work/queue"
  fi
  if [ "$name" != "$realname" ] && [ ! -e "$fw/$name" ]; then
    ln -s "$realname" "$fw/$name"
  fi
}

for f in "$exe_dir"/*; do
  if [ -f "$f" ] && [ ! -L "$f" ] && is_macho "$f"; then
    printf '%s\t%s\n' "$f" "$exe_dir" >>"$work/queue"
  fi
done

i=1
while :; do
  line=$(sed -n "${i}p" "$work/queue")
  [ -n "$line" ] || break
  i=$((i + 1))
  f=${line%%$'\t'*}
  src=${line#*$'\t'}
  for r in $(refs "$f"); do
    if is_system "$r"; then continue; fi
    name=$(basename "$r")
    path=$(resolve "$r" "$f" "$src")
    if [ -z "$path" ]; then
      echo "error: $(basename "$f") loads $r, which is not in $prefix" >&2
      exit 1
    fi
    case "$r" in
      @rpath/*.framework/*)
        bundle_framework "$path" "${r#@rpath/}"
        continue
        ;;
    esac
    bundle "$path" "$name"
    # Everything in Frameworks is found through @rpath, whatever form the
    # reference had in the prefix.
    if [ "$r" != "@rpath/$name" ]; then
      "$INSTALL_NAME_TOOL" -change "$r" "@rpath/$name" "$f"
      echo "$f" >>"$work/modified"
    fi
  done
done

# The executables look in Frameworks, and nowhere absolute: an LC_RPATH into
# the build machine's prefix would let a copy of the library there win over
# the bundled one on that machine, and hide a missing library.
for f in "$exe_dir"/*; do
  if [ ! -f "$f" ] || [ -L "$f" ] || ! is_macho "$f"; then continue; fi
  have=no
  for p in $(rpaths "$f"); do
    case "$p" in
      @executable_path/../Frameworks) have=yes ;;
      /*)
        "$INSTALL_NAME_TOOL" -delete_rpath "$p" "$f"
        echo "$f" >>"$work/modified"
        ;;
    esac
  done
  if [ "$have" = no ]; then
    "$INSTALL_NAME_TOOL" -add_rpath @executable_path/../Frameworks "$f"
    echo "$f" >>"$work/modified"
  fi
done

sort -u "$work/modified" | while read -r f; do
  "$CODESIGN" --force --sign - "$f"
done

n=$(grep -c -v "^$exe_dir/" "$work/queue" || true)
size=$(du -sh "$fw" | cut -f1)
echo "bundled $n libraries ($size) into $fw"
