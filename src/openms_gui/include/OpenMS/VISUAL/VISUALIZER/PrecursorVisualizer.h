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
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  /**
      @brief Class that displays all meta information for Precursor objects

      This class provides all functionality to view the meta information of an object of type Precursor.
  */
  class OPENMS_GUI_DLLAPI PrecursorVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<Precursor>
  {
    Q_OBJECT

public:

    ///Constructor
    PrecursorVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:
    ///@name Edit fields and buttons
    //@{
    QLineEdit * mz_;
    QLineEdit * int_;
    QLineEdit * charge_;
    QLineEdit * window_up_;
    QLineEdit * window_low_;
    QListWidget * activation_methods_;
    QLineEdit * activation_energy_;
    //@}

    //Docu in base class
    void update_() override;
  };
}
