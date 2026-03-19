#!/usr/bin/env bash
set -euo pipefail

yum install -y openssh-clients

yum install -y --setopt=install_weak_deps=False \
  libcurl-devel \
  dpkg \
  dpkg-dev

is_pkg_available() {
  local pkg="$1"
  yum -q list available "$pkg" >/dev/null 2>&1 || yum -q list installed "$pkg" >/dev/null 2>&1
}

install_first_available() {
  local group_name="$1"
  shift
  local pkg

  for pkg in "$@"
  do
    if is_pkg_available "$pkg"
    then
      yum install -y --setopt=install_weak_deps=False "$pkg"
      return 0
    fi
  done

  echo "Skipping ${group_name}: no matching package available (${*})."
  return 0
}

install_first_available "Arrow devel" \
  arrow-devel \
  libarrow-devel

install_first_available "Arrow compute devel" \
  arrow-compute-devel \
  libarrow-compute-devel \
  libarrow-devel

install_first_available "Parquet devel" \
  parquet-devel \
  libparquet-devel \
  parquet-libs-devel
