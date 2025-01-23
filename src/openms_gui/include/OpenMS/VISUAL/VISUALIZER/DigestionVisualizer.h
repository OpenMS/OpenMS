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
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>
#include <OpenMS/METADATA/Digestion.h>


namespace OpenMS
{
  /**
      @brief Class that displays all meta information of digestion objects.

      This class provides all functionality to view the meta information of an object of type Digestion.
  */
  class OPENMS_GUI_DLLAPI DigestionVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<Digestion>
  {
    Q_OBJECT

public:

    ///Constructor
    DigestionVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * treatmenttype_;
    QTextEdit * treatmentcomment_;
    QLineEdit * digestionenzyme_;
    QLineEdit * digestiontime_;
    QLineEdit * digestiontemperature_;
    QLineEdit * digestionPH_;
    //@}

    //Docu in base class
    void update_() override;
  };

}
