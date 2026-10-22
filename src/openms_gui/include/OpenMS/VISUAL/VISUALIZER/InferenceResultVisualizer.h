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
  class MetaDataBrowser;

  /**
      @brief Class that displays the proteins of an inference result of identification data (its protein run)

      The score type, its orientation and the significance threshold of the proteins can be edited; the runs that
      the inference used are shown. Its proteins can be filtered by score in the tree of the browser.
  */
  class OPENMS_GUI_DLLAPI InferenceResultVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<ProteinIdentification>
  {
    Q_OBJECT

public:

    ///Constructor
    InferenceResultVisualizer(bool editable = false, QWidget * parent = nullptr, MetaDataBrowser * caller = nullptr);

    /// Loads the proteins of @p result (a copy that the browser edits), which has the item @p tree_item_id in the browser
    void load(IdentificationData::InferenceResult & result, int tree_item_id);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_() override;

    /// Shows the proteins with a score at least as good as the threshold
    void updateTree_();

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * identifier_;
    QTextEdit * inputs_;
    QLineEdit * score_type_;
    QComboBox * higher_better_;
    QLineEdit * significance_threshold_;
    QLineEdit * filter_threshold_;
    //@}

    //Docu in base class
    void update_() override;

private:
    /// The browser that shows the result
    MetaDataBrowser * caller_;
    /// The item of the result in the tree of the browser
    int tree_id_ = 0;
  };
}
