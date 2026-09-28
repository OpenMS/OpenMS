Libraries from the OpenMS contrib
=================================

A build against the OpenMS contrib (https://github.com/OpenMS/contrib), such as the pyOpenMS
wheels, links the libraries below from it. The contrib builds them from the source archives at
https://github.com/OpenMS/contrib-sources/releases/tag/3.6.0. Each folder next to this file holds
the license and notice texts of one library, as the library's own source distribution has them;
a build installs the folders of the libraries it took from the contrib.

Library     Version  License                                 Source code
----------  -------  --------------------------------------  ----------------------------------------------------------
boost       1.87.0   BSL-1.0                                 https://github.com/boostorg/boost/tree/boost-1.87.0
bzip2       1.0.5    bzip2-1.0.6                             https://sourceware.org/git/bzip2.git, tag bzip2-1.0.5
zlib        1.3.1    Zlib                                    https://github.com/madler/zlib/tree/v1.3.1
curl        8.12.1   curl                                    https://github.com/curl/curl/tree/curl-8_12_1
eigen       3.4.0    MPL-2.0, parts Apache-2.0, BSD-3-Clause, https://gitlab.com/libeigen/eigen/-/tree/3.4.0
                     Minpack (see COPYING.README)
libsvm      3.12     BSD-3-Clause                            https://github.com/cjlin1/libsvm/tree/v312
libzip      1.11.4   BSD-3-Clause                            https://github.com/nih-at/libzip/tree/v1.11.4
xerces-c    3.2.0    Apache-2.0                              https://github.com/apache/xerces-c/tree/v3.2.0
coinmp      1.8.3    Cbc 2.9.6, Cgl 0.59.7, Clp 1.16.8,      https://www.coin-or.org/download/source/CoinMP/CoinMP-1.8.3.tgz
                     CoinUtils 2.10.10, Osi 0.107.6: EPL-1.0
arrow       23.0.0   Apache-2.0                              https://github.com/apache/arrow/tree/apache-arrow-23.0.0

Arrow builds and links these itself:

snappy      1.2.2    BSD-3-Clause                            https://github.com/google/snappy/tree/1.2.2
zstd        1.5.7    BSD-3-Clause OR GPL-2.0-only            https://github.com/facebook/zstd/tree/v1.5.7
thrift      0.22.0   Apache-2.0                              https://github.com/apache/thrift/tree/v0.22.0
xsimd       14.0.0   BSD-3-Clause                            https://github.com/xtensor-stack/xsimd/tree/14.0.0
rapidjson   232389d  MIT                                     https://github.com/Tencent/rapidjson/tree/232389d4f1012dddec4ef84861face2d2ba85709
mimalloc    3.1.5    MIT                                     https://github.com/microsoft/mimalloc/tree/v3.1.5

The source code of Eigen (MPL-2.0) and of the COIN-OR projects (EPL-1.0) is available at the
addresses above; the COIN-OR projects are also at https://github.com/coin-or/<project>, tag
releases/<version>. zstd is used under its BSD-3-Clause license. Of the COIN-OR libraries,
OpenMS links Cbc, Cgl, Clp, CoinUtils and Osi, not CoinMP itself.

Eigen's NonLinearOptimization module, which OpenMS uses, is a port of MINPACK:
This product includes software developed by the University of Chicago, as Operator of Argonne
National Laboratory.
