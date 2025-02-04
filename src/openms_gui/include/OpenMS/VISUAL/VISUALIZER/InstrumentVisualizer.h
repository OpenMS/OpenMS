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
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  /**
      @brief Class that displays all meta information for an MS instrument

      This class provides all functionality to view the meta information of an object of type Instrument.
  */
  class OPENMS_GUI_DLLAPI InstrumentVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<Instrument>
  {
    Q_OBJECT

public:

    ///Constructor
    InstrumentVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * name_;
    QLineEdit * vendor_;
    QLineEdit * model_;
    QTextEdit * customizations_;
    QComboBox * ion_optics_;
    //@}

    //Docu in base class
    void update_() override;
  };
}
