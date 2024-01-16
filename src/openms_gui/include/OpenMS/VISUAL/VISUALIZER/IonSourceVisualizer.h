// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  /**
      @brief Class that displays all meta information for IonSource objects

      This class provides all functionality to view the meta information of an object of type IonSource
  */
  class OPENMS_GUI_DLLAPI IonSourceVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<IonSource>
  {
    Q_OBJECT

public:

    ///Constructor
    IonSourceVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * order_;
    QComboBox * inlet_type_;
    QComboBox * ionization_method_;
    QComboBox * polarity_;
    //@}

    //Docu in base class
    void update_() override;
  };
}
