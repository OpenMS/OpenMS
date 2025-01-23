// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

//QT

namespace OpenMS
{
  /**
      @brief Class that displays all meta information for AcquisitionInfo objects

      This class provides all functionality to view the meta information of an object of type AcquisitionInfo.
  */
  class OPENMS_GUI_DLLAPI AcquisitionInfoVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<AcquisitionInfo>
  {
    Q_OBJECT

public:

    ///Constructor
    AcquisitionInfoVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    /// Edit field for the method
    QLineEdit * acquisitioninfo_method_;

    //Docu in base class
    void update_() override;
  };

}
