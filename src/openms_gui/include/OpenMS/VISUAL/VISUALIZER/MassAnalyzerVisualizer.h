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
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  /**
      @brief Class that displays all meta information for MassAnalyzer objects

      This class provides all functionality to view the meta information of an object of type MassAnalyzer.
  */
  class OPENMS_GUI_DLLAPI MassAnalyzerVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<MassAnalyzer>
  {
    Q_OBJECT

public:

    ///Constructor
    MassAnalyzerVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    ///@name edit fields to modify properties
    //@{
    QLineEdit * order_;
    QLineEdit * res_;
    QLineEdit * acc_;
    QLineEdit * scan_rate_;
    QLineEdit * scan_time_;
    QLineEdit * TOF_;
    QLineEdit * iso_;
    QLineEdit * final_MS_;
    QLineEdit * magnetic_fs_;
    QComboBox * type_;
    QComboBox * res_method_;
    QComboBox * res_type_;
    QComboBox * scan_dir_;
    QComboBox * scan_law_;
    QComboBox * reflectron_state_;
    //@}

    //Docu in base class
    void update_() override;
  };
}
