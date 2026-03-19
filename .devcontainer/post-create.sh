#!/usr/bin/env bash
set -euo pipefail

yum install -y openssh-clients

yum install -y --setopt=install_weak_deps=False \
  libcurl-devel \
  dpkg \
  dpkg-dev

if ! rpm -q apache-arrow-release >/dev/null 2>&1
then
  yum install -y https://packages.apache.org/artifactory/arrow/almalinux/9/apache-arrow-release-latest.rpm
fi

# EPEL Arrow/Parquet devel packages conflict with Apache Arrow 23 package layout.
# Remove them first so Arrow 23 installation can proceed cleanly.
yum remove -y libarrow-devel parquet-libs-devel libparquet-devel 2>/dev/null || true

yum install -y --setopt=install_weak_deps=False --allowerasing \
  arrow-devel-23* \
  arrow-compute-devel-23* \
  parquet-devel-23*
