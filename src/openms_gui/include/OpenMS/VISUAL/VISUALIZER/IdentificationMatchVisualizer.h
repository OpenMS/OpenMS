// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  /**
      @brief Class that displays a match (candidate) of an identification of identification data

      The molecule, charge, target/decoy state, scores and sequence evidence of the match are shown (metadata is
      edited in the metadata view of the match).
  */
  class OPENMS_GUI_DLLAPI IdentificationMatchVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<IdentificationData::MatchData>
  {
    Q_OBJECT

public:

    ///Constructor
    IdentificationMatchVisualizer(bool editable = false, QWidget * parent = nullptr);

    /// Loads the data of @p match of @p run (with the scores that the run holds for it)
    void load(IdentificationData::MatchData & data, const IdentificationData::Run & run, const IdentificationData::Match & match);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_() override;

protected:

    ///@name Edit fields and buttons
    //@{
    QTextEdit * molecule_;
    QLineEdit * charge_;
    QLineEdit * target_decoy_;
    QTextEdit * scores_;
    QTextEdit * evidence_;
    //@}

    //Docu in base class
    void update_() override;
  };
}
