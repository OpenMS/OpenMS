// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $

int pch_probe();
extern "C" int pch_c_probe(void);

int main()
{ return pch_probe() == pch_c_probe() ? 0 : 1; }
