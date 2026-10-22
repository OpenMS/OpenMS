// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <QtTest/QtTest>
#include <QtGui>

namespace OpenMS
{
  /// The identification data of maps in the MetaDataBrowser: shown, edited and written back
  class TestMetaDataBrowser : public QObject
  {
    Q_OBJECT

    private slots:
      void testIdentificationTree();

      void testFeatureIdentifications();

      void testEditAndSave();
  };
}
