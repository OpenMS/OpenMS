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
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

class QListWidget;

namespace OpenMS
{
  /**
      @brief Class that displays all meta information for DataProcessing objects

      This class provides all functionality to view the meta information of an object of type DataProcessing.
  */
  class OPENMS_GUI_DLLAPI DataProcessingVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<DataProcessing>
  {
    Q_OBJECT

public:

    /// Constructor
    DataProcessingVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * completion_time_;
    QListWidget * actions_;
    //@}

    //Docu in base class
    void update_() override;
  };
}
