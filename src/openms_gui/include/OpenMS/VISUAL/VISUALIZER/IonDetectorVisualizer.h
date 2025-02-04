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
#include <OpenMS/METADATA/IonDetector.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  /**
      @brief Class that displays all meta information for IonDetector objects

      This class provides all functionality to view the meta information of an object of type IonDetector.
  */
  class OPENMS_GUI_DLLAPI IonDetectorVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<IonDetector>
  {
    Q_OBJECT

public:

    ///Constructor
    IonDetectorVisualizer(bool editable = false, QWidget * parent = nullptr);

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
    QLineEdit * freq_;
    QComboBox * type_;
    QComboBox * ac_mode_;
    //@}

    //Docu in base class
    void update_() override;
  };
}
