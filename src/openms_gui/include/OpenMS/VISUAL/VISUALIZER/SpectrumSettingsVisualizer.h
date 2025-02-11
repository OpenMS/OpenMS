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

//OpenMS
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

class QTextEdit;

namespace OpenMS
{
  /**
      @brief Class that displays all meta information for SpectrumSettings objects

      This class provides all functionality to view the meta information of an object of type SpectrumSettings.
  */
  class OPENMS_GUI_DLLAPI SpectrumSettingsVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<SpectrumSettings>
  {
    Q_OBJECT

public:

    ///Constructor
    SpectrumSettingsVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:
    /// The date of this experiment
    QLineEdit * native_id_;
    /// The type of this experiment
    QComboBox * type_;
    /// The date of this experiment
    QTextEdit * comment_;

    //Docu in base class
    void update_() override;
  };
}
