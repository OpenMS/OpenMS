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

#include <OpenMS/METADATA/Modification.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  /**
      @brief Class that displays all meta information of modification objects.

      This class provides all functionality to view the meta information of an object of type Modification.
  */
  class OPENMS_GUI_DLLAPI ModificationVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<Modification>
  {
    Q_OBJECT

public:

    ///Constructor
    ModificationVisualizer(bool editable = false, QWidget * parent = nullptr);

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
    QLineEdit * modificationname_;
    QLineEdit * modificationmass_;
    QComboBox * modificationspecificity_;
    QLineEdit * modificationAA_;
    //@}

    //Docu in base class
    void update_() override;
  };

}
