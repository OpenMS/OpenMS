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
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  /**
      @brief Class that displays all meta information for ContactPerson objects

      This class provides all functionality to view the meta information of an object of type ContactPerson.
  */
  class OPENMS_GUI_DLLAPI ContactPersonVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<ContactPerson>
  {
    Q_OBJECT

public:

    ///Constructor
    ContactPersonVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * firstname_;
    QLineEdit * lastname_;
    QLineEdit * institution_;
    QLineEdit * email_;
    QLineEdit * contact_info_;
    QLineEdit * address_;
    QLineEdit * url_;
    //@}

    //Docu in base class
    void update_() override;
  };

}
