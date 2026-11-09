#!/usr/bin/env python3
"""bundle-licences.py PREFIX DEST --libs DIR [--libs DIR ...] [--require]

Ship the licence texts of the third-party libraries a package bundles.

F2Dock's packages carry shared libraries from the cvcpkg prefix they were
built against: the Linux tarball's lib/, the Windows zip's bin/ (every DLL in
the prefix), the macOS .app's Contents/Frameworks. Their licences require the
texts to travel with them -- gmp and mpfr are LGPL-3.0, libstdc++ is GPL-3.0
with the GCC Runtime Library Exception, most of the rest want their copyright
notice reproduced. cvcpkg puts each package's upstream texts in
<prefix>/share/licenses/<package>/.

This copies every <prefix>/share/licenses/<package>/ into DEST/<package>/,
writes DEST/README.txt listing each bundled library's package, version and
licence, and reports every bundled library whose package has no licence
texts. With --require that is an error: a package with a library it cannot
license is not one to publish.

A library's package is found through the prefix's installed manifests
(<prefix>/share/libcvc-deps/<package>/manifest.yaml), which list every file a
package installed. A library no manifest lists -- its bundle was published
without one -- is reported as such.

Runs on the Python that installed cvcpkg, so PyYAML is there.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import yaml

# What counts as a shared library, by platform spelling. Symlinks are names
# of a library, not separate libraries; they are matched by name too.
SUFFIXES = (".dylib", ".dll", ".so")


def is_library(name: str) -> bool:
    return name.endswith(SUFFIXES) or ".so." in name or (".framework" in name)


def load_manifests(prefix: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    root = prefix / "share" / "libcvc-deps"
    if not root.is_dir():
        return out
    for mf in sorted(root.glob("*/manifest.yaml")):
        with open(mf, encoding="utf-8") as fh:
            out[mf.parent.name] = yaml.safe_load(fh) or {}
    return out


def owners_by_basename(manifests: dict[str, dict]) -> dict[str, set[str]]:
    owners: dict[str, set[str]] = {}
    for pkg, m in manifests.items():
        for f in (m.get("contents") or {}).get("files") or []:
            name = f.rstrip("/").split("/")[-1]
            owners.setdefault(name, set()).add(pkg)
            # A framework is listed file by file; its top directory names it.
            for part in f.split("/"):
                if part.endswith(".framework"):
                    owners.setdefault(part, set()).add(pkg)
    return owners


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("prefix", type=Path)
    ap.add_argument("dest", type=Path)
    ap.add_argument("--libs", type=Path, action="append", default=[], required=True)
    ap.add_argument("--require", action="store_true")
    args = ap.parse_args()

    prefix = args.prefix.resolve()
    lic_root = prefix / "share" / "licenses"
    manifests = load_manifests(prefix)
    owners = owners_by_basename(manifests)

    args.dest.mkdir(parents=True, exist_ok=True)
    copied: set[str] = set()
    if lic_root.is_dir():
        for d in sorted(p for p in lic_root.iterdir() if p.is_dir()):
            if any(d.iterdir()):
                shutil.copytree(d, args.dest / d.name, dirs_exist_ok=True)
                copied.add(d.name)

    bundled: dict[str, set[str]] = {}  # package -> library names
    unowned: list[str] = []
    for libdir in args.libs:
        for p in sorted(libdir.iterdir()):
            if not is_library(p.name):
                continue
            pkgs = owners.get(p.name)
            if not pkgs:
                unowned.append(p.name)
                continue
            for pkg in pkgs:
                bundled.setdefault(pkg, set()).add(p.name)

    missing = sorted(pkg for pkg in bundled if pkg not in copied)

    lines = [
        "Licence texts of the third-party libraries bundled with this package.",
        "",
        "Each directory here holds the upstream licence texts of one cvcpkg",
        "package, as that package installed them in share/licenses/<package>/.",
        "",
        f"{'package':<26} {'version':<24} licence",
        "-" * 78,
    ]
    for pkg in sorted(bundled):
        meta = manifests.get(pkg, {})
        version = (meta.get("bundle") or {}).get("version", "?")
        licence = (meta.get("meta") or {}).get("license", "?")
        mark = "" if pkg in copied else "   [NO LICENCE TEXT IN THE BUNDLE]"
        lines.append(f"{pkg:<26} {version:<24} {licence}{mark}")
    if unowned:
        lines += ["", "Bundled libraries no installed manifest lists:"]
        lines += [f"  {n}" for n in unowned]
    (args.dest / "README.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(
        f"licences: {len(bundled)} packages bundled, {len(copied)} licence directories "
        f"copied to {args.dest}"
    )
    for pkg in missing:
        libs = ", ".join(sorted(bundled[pkg])[:3])
        print(f"{'ERROR' if args.require else 'WARNING'}: {pkg} has no licence texts ({libs})")
    if unowned:
        print(
            f"{'ERROR' if args.require else 'WARNING'}: {len(unowned)} bundled libraries "
            f"have no owning manifest, e.g. {', '.join(unowned[:3])}"
        )
    if args.require and (missing or unowned):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
