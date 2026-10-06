// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#pragma once

namespace SQLite
{
class Database;
}
namespace OpenMS
{
class IdentificationData;
}

namespace OpenMS::Internal
{
/// Direct SQLite persistence for the owning identification model (OMS schema 6).
void storeOMSIdentifications(SQLite::Database& db, const IdentificationData& data);
void loadOMSIdentifications(SQLite::Database& db, IdentificationData& data);
} // namespace OpenMS::Internal
