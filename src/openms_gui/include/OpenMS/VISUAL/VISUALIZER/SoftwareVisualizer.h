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
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{

  class MetaDataBrowser;

  /**
      @brief Class that displays all meta information for Software objects

      This class provides all functionality to view the meta information of an object of type Software.
  */
  class OPENMS_GUI_DLLAPI SoftwareVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<Software>
  {
    Q_OBJECT

public:

    ///Constructor
    SoftwareVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * software_name_;
    QLineEdit * software_version_;
    //@}

    //Docu in base class
    void update_() override;
  };
}
