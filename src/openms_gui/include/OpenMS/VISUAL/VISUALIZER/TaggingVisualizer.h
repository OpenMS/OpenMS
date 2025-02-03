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
#include <OpenMS/METADATA/Tagging.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

class QDoubleValidator;

namespace OpenMS
{
  /**
      @brief Class that displays all meta information of tagging objects.

      This class provides all functionality to view the meta information of an object of type Tagging.
  */
  class OPENMS_GUI_DLLAPI TaggingVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<Tagging>
  {
    Q_OBJECT

public:

    ///Constructor
    TaggingVisualizer(bool editable = false, QWidget * parent = nullptr);

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
    QLineEdit * taggingmass_shift_;
    QComboBox * taggingvariant_;
    //@}

    //Docu in base class
    void update_() override;
  };

}
