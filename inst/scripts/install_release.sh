#!/usr/bin/env bash
# Installs an earlier prolfqua release and the prolfquapp of its time into results/releases/<name>/lib, for
# `make score-history-release`; the other packages come from the normal library. The source trees go to
# results/releases/<name>/src/<package>-<ref>, the install logs next to them.
#   inst/scripts/install_release.sh <name> <prolfqua commit or tag> <prolfquapp commit or tag>
# Run from the prolfquabenchmark root, next to ../prolfqua and ../prolfquapp.
set -euo pipefail

dir=$(pwd)/results/releases/$1
mkdir -p "$dir/lib"
install_package() {
  local package=$1 ref=$2
  local src=$dir/src/$package-$ref
  mkdir -p "$src"
  git -C "../$package" archive "$ref" | tar -x -C "$src"
  R_LIBS=$dir/lib R CMD INSTALL -l "$dir/lib" "$src" > "$dir/install_$package.log" 2>&1
  echo "$package $ref installed into $dir/lib"
}
install_package prolfqua "$2"
install_package prolfquapp "$3"
