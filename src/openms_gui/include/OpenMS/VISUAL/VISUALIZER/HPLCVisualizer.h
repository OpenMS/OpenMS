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
#include <OpenMS/METADATA/HPLC.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>


namespace OpenMS
{
  /**
      @brief Class that displays all meta information for HPLC objects

      This class provides all functionality to view the meta information of an object of type HPLC.
  */
  class OPENMS_GUI_DLLAPI HPLCVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<HPLC>
  {
    Q_OBJECT

public:

    ///Constructor
    HPLCVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * hplcinstrument_;
    QLineEdit * hplccolumn_;
    QLineEdit * hplctemperature_;
    QLineEdit * hplcpressure_;
    QLineEdit * hplcflux_;
    QTextEdit * hplccomment_;
    //@}

    //Docu in base class
    void update_() override;
  };

}
