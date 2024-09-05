// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  /**
      @brief Class that displays all meta information for PeptideHit objects

      This class provides all functionality to view the meta information of an object of type PeptideHit.
  */
  class OPENMS_GUI_DLLAPI PeptideHitVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<PeptideHit>
  {
    Q_OBJECT

public:

    ///Constructor
    PeptideHitVisualizer(bool editable = false, QWidget * parent = nullptr);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * peptidehit_score_;
    QLineEdit * peptidehit_charge_;
    QLineEdit * peptidehit_rank_;
    QTextEdit * peptidehit_sequence_;
    //@}

    //Docu in base class
    void update_() override;
  };


}
