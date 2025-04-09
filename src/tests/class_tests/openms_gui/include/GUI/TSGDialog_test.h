// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <QtTest/QtTest>
#include <QtGui>
#include <QtWidgets/QSpinBox>

#include <OpenMS/VISUAL/DIALOGS/TheoreticalSpectrumGenerationDialog.h>
#include <ui_TheoreticalSpectrumGenerationDialog.h>

#define UI dialog_.ui_

namespace OpenMS
{
  class TestTSGDialog : public QObject
  {
    Q_OBJECT

    public:
    TestTSGDialog() : dialog_() {}

    ~TestTSGDialog()
    {
      dialog_.destroy();
    }

    private slots:
      void testConstruction();
      
      void testGui();

      void testParameterImport();

      void testSpectrumCalculation();

      void testErrors();

    private:
      template<typename T> // template for QSpinBox and QDoubleSpinBox
      void testSpinBox_(T* box, std::string str_value = "2");

      void testIonsIntensities_();

      void testSequenceInput_(QString input);

      void testIsotopeModel_(bool skip_none = false);

      void checkMessageBoxExists_();

      void testMessageBoxes_();

      TheoreticalSpectrumGenerationDialog dialog_;
  };
}