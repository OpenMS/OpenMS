// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{

  class MetaDataBrowser;

  /**
      @brief Class that displays all meta information for ExperimentalSettings objects

      This class provides all functionality to view the meta information of an object of type ExperimentalSettings.
  */
  class OPENMS_GUI_DLLAPI ExperimentalSettingsVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<ExperimentalSettings>
  {
    Q_OBJECT

public:

    ///Constructor
    ExperimentalSettingsVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:
    /// The date of this experiment
    QLineEdit * datetime_;
    /// The comment to this experiment
    QTextEdit * comment_;
    /// Fraction identifier
    QLineEdit * fraction_identifier_;

    //Docu in base class
    void update_() override;
  };
}
