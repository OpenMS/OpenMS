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
      @brief Class that displays an identification of identification data (its observation, e.g. of a spectrum)

      The data reference, RT and m/z of the observation and the primary score of its run are shown (metadata is edited
      in the metadata view of the identification). Its matches can be filtered by score in the tree of the browser.
  */
  class OPENMS_GUI_DLLAPI IdentificationVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<IdentificationData::Observation>
  {
    Q_OBJECT

public:

    ///Constructor
    IdentificationVisualizer(bool editable = false, QWidget * parent = nullptr, MetaDataBrowser * caller = nullptr);

    /// Loads the observation of an identification of @p run, which has the item @p tree_item_id in the browser
    void load(IdentificationData::Observation & observation, const IdentificationData::Run & run, int tree_item_id);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_() override;

    /// Shows the matches with a score at least as good as the threshold
    void updateTree_();

protected:

    ///@name Edit fields and buttons
    //@{
    QLineEdit * run_;
    QLineEdit * data_id_;
    QLineEdit * rt_;
    QLineEdit * mz_;
    QLineEdit * score_;
    QLineEdit * filter_threshold_;
    //@}

    //Docu in base class
    void update_() override;

private:
    /// The browser that shows the identification
    MetaDataBrowser * caller_;
    /// The item of the identification in the tree of the browser
    int tree_id_ = 0;
    /// Whether higher primary scores are better
    bool higher_better_ = true;
  };
}
