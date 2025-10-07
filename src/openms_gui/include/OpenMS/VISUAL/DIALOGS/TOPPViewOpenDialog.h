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

#include <QtWidgets/QDialog>

#include <map>
class QAbstractButton;

namespace Ui
{
  class TOPPViewOpenDialogTemplate;
}

namespace OpenMS
{
  class Param;
  class String;
  /**
      @brief Dataset opening options for TOPPView

      @ingroup TOPPView_elements
  */
  class OPENMS_GUI_DLLAPI TOPPViewOpenDialog :
    public QDialog
  {
    Q_OBJECT

public:
    /// Constructor
    TOPPViewOpenDialog(const String & data_name, bool as_window, bool as_2d, bool cutoff, QWidget * parent = nullptr);
    /// Destructor
    ~TOPPViewOpenDialog() override;

    /// Returns true, if 2D mode is to be used for maps
    bool viewMapAs2D() const;
    /// Returns true, if 1D mode is to be used for maps
    bool viewMapAs1D() const;
    /// Returns if the low intensity peaks should be hidden
    bool isCutoffEnabled() const;
    /// Returns if the data is DIA / SWATH-MS data
    bool isDataDIA() const;
    /// Returns true, if the data should be opened in a new window
    bool openAsNewWindow() const;
    ///Returns the index of the selected merge layer. If the option is not selected -1 is returned.
    Int getMergeLayer() const;

    /// Disables view dimension section and sets the selected option
    void disableDimension(bool as_2d);
    /// Disables cutoff section and sets the selected option
    void disableCutoff(bool cutoff_on);
    /// Disables opening location section and sets the selected option
    void disableLocation(bool window);
    /**
        @brief Sets the possible merge layers (index and name) and activates the option

        It is deactivated by default and can be deactivated manually by passing an empty list.
    */
    void setMergeLayers(const std::map<Size, String> & layers);

protected slots:
    ///slot that disables 2D/3D options, when as layer is selected
    void updateViewMode_(QAbstractButton * button);

protected:
    ///Stores if this option is disabled, to avoid activating it in updateViewMode_()
    bool map_as_2d_disabled_;

private: 
    Ui::TOPPViewOpenDialogTemplate* ui_;
  };

}
