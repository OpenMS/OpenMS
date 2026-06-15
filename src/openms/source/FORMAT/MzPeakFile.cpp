// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/MzPeakFile.h>

namespace OpenMS
{

MzPeakFile::MzPeakFile() = default;

MzPeakFile::~MzPeakFile() = default;

void MzPeakFile::load(const String& /* filename */, MapType& /* map */) const
{ throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }

void MzPeakFile::store(const String& /* filename */, const MapType& /* map */) const
{ throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }

void MzPeakFile::transform(const String& /* filename_in */,
                           Interfaces::IMSDataConsumer* /* consumer */,
                           bool /* skip_full_count */,
                           bool /* skip_first_pass */) const
{ throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }

} // namespace OpenMS
