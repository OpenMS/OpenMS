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
      @brief Class that displays the settings of a run of identification data (IdentificationData::RunSettings)

      The software, its version, the date and the search parameters can be edited. The identifier, molecule kind,
      databases and score definitions of the run are shown (see setRun()).
  */
  class OPENMS_GUI_DLLAPI IdentificationRunVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<IdentificationData::RunSettings>
  {
    Q_OBJECT

public:

    ///Constructor
    IdentificationRunVisualizer(bool editable = false, QWidget * parent = nullptr);

    /// Shows what the run has besides its settings (identifier, molecule kind, databases, score definitions)
    void setRun(const IdentificationData::Run& run);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_() override;

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * identifier_;
    QLineEdit * kind_;
    QTextEdit * databases_;
    QTextEdit * scores_;
    QLineEdit * software_;
    QLineEdit * software_version_;
    QLineEdit * date_;
    QLineEdit * charges_;
    QLineEdit * missed_cleavages_;
    QLineEdit * fragment_tolerance_;
    QComboBox * fragment_tolerance_ppm_;
    QLineEdit * precursor_tolerance_;
    QComboBox * precursor_tolerance_ppm_;
    QComboBox * mass_type_;
    QLineEdit * enzyme_;
    //@}

    //Docu in base class
    void update_() override;
  };
}
