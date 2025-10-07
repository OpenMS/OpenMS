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

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/VISUAL/PlotCanvas.h> // for AreaXYType

#include <QtWidgets/QDialog>

namespace Ui
{
  class Plot2DGoToDialogTemplate;
}

namespace OpenMS
{

  class String;

  /**
      @brief GoTo dialog used to zoom to a m/z and retention time range or to a feature.

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI Plot2DGoToDialog :
    public QDialog
  {
    Q_OBJECT

public:
    using AreaXYType = PlotCanvas::AreaXYType;

    ///Constructor
    /// @param parent Parent widget
    /// @param x_name Name of the x_axis dimension
    /// @param y_name Name of the y_axis dimension
    Plot2DGoToDialog(QWidget* parent, std::string_view x_name, std::string_view y_name);
    ///Destructor
    ~Plot2DGoToDialog() override;

    /// Returns if a feature UID was set an a feature should be displayed (false), otherwise, show a range (true)
    bool showRange() const;

    bool checked();

    ///@name Methods for ranges
    //@{
    ///Sets the data range to display initially
    void setRange(const AreaXYType& range);
    ///Sets the data range of the complete experiment for better navigation with the dialog
    void setMinMaxOfRange(const AreaXYType& max_range);

    /// Query the range set by the user.
    /// If any dimension is <1, it is extended to at least 1 to ensure proper displaying.
    AreaXYType getRange();
    //@}

    ///@name Methods for feature numbers
    //@{
    ///Returns the selected feature numbers. If a number is returned, the feature rather than the range should be displayed.
    String getFeatureNumber() const;
    ///Disables the feature number field
    void enableFeatureNumber(bool);
    //@}

  private:
    Ui::Plot2DGoToDialogTemplate* ui_;

  };

}
