#!/usr/bin/env bash
set -euo pipefail

yum install -y openssh-clients

if rpm -q libarrow-devel >/dev/null 2>&1 || rpm -q libarrow >/dev/null 2>&1
then
  yum remove -y libarrow-devel libarrow || true
fi

yum install -y --setopt=install_weak_deps=False \
  libcurl-devel \
  arrow-devel \
  arrow-compute-devel \
  parquet-devel \
  dpkg \
  dpkg-dev
