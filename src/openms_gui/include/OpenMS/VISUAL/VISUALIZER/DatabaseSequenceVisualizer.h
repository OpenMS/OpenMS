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
      @brief Class that displays a database sequence of a run of identification data (e.g. a protein)

      The accession, database, target/decoy state, description and sequence are shown (metadata is edited in the
      metadata view of the sequence).
  */
  class OPENMS_GUI_DLLAPI DatabaseSequenceVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<IdentificationData::DatabaseSequence>
  {
    Q_OBJECT

public:

    ///Constructor
    DatabaseSequenceVisualizer(bool editable = false, QWidget * parent = nullptr);

    /// Loads a database sequence of @p run
    void load(IdentificationData::DatabaseSequence & sequence, const IdentificationData::Run & run);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_() override;

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * accession_;
    QLineEdit * database_;
    QLineEdit * target_decoy_;
    QTextEdit * description_;
    QTextEdit * sequence_;
    //@}

    //Docu in base class
    void update_() override;
  };
}
