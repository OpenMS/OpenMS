// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QtWidgets/QDialog>

#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/VISUAL/HistogramWidget.h>

namespace OpenMS
{
  /**
      @brief Dialog that show a HistogramWidget.

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI HistogramDialog :
    public QDialog
  {
    Q_OBJECT

public:
    /// Constructor
    HistogramDialog(const Math::Histogram<> & distribution, QWidget * parent = nullptr);
    /// Destructor
    ~HistogramDialog() override;

    /// Returns the value of the left splitter
    float getLeftSplitter();
    /// Returns the value of the right splitter
    float getRightSplitter();

    /// Sets the value of the left splitter
    void setLeftSplitter(float position);
    /// Sets the value of the right splitter
    void setRightSplitter(float position);

    /// Sets the axis legend
    void setLegend(const String & legend);
    /// Sets log mode
    void setLogMode(bool log_mode);

protected:
    HistogramWidget * mw_;
  };

} //namespace

