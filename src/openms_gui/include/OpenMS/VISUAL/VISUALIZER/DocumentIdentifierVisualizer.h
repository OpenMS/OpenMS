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
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

//QT

namespace OpenMS
{
  /**
      @brief Class that displays all meta information for DocumentIdentifier objects

      This class provides all functionality to view the meta information of an object of type DocumentIdentifier.
  */
  class OPENMS_GUI_DLLAPI DocumentIdentifierVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<DocumentIdentifier>
  {
    Q_OBJECT

public:

    ///Constructor
    DocumentIdentifierVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:
    ///@name Edit fields and buttons
    //@{
    QLineEdit * identifier_;
    QLineEdit * file_path_;
    QLineEdit * file_type_;
    //@}

    //Docu in base class
    void update_() override;
  };
}
