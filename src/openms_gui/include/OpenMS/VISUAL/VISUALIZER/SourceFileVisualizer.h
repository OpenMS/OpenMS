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
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  /**
      @brief Class that displays all meta information for SourceFile objects

      This class provides all functionality to view the meta information of an object of type SourceFile.
  */
  class OPENMS_GUI_DLLAPI SourceFileVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<SourceFile>
  {
    Q_OBJECT

public:

    ///Constructor
    SourceFileVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * name_of_file_;
    QLineEdit * path_to_file_;
    QLineEdit * file_size_;
    QLineEdit * file_type_;
    QLineEdit * checksum_;
    QComboBox * checksum_type_;
    QLineEdit * native_id_type_;
    //@}

    //Docu in base class
    void update_() override;
  };

}
