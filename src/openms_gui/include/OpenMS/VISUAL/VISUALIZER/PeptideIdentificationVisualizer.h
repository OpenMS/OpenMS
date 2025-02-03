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
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  class MetaDataBrowser;

  /**
      @brief Class that displays all meta information for PeptideIdentification objects

      This class provides all functionality to view the meta information of an object of type PeptideIdentification.
  */
  class OPENMS_GUI_DLLAPI PeptideIdentificationVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<PeptideIdentification>
  {
    Q_OBJECT

public:
    ///Constructor
    PeptideIdentificationVisualizer(bool editable = false, QWidget * parent = nullptr, MetaDataBrowser * caller = nullptr);

    /// Loads the meta data from the object to the viewer. Gets the id of the item in the tree as parameter.
    void load(PeptideIdentification & s, int tree_item_id);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    ///Undo the changes made in the GUI.
    void undo_();

    /**
        @brief Updates the tree by calling MetaDataBrowser::updatePeptideHits(PeptideIdentification, int)

        Calls MetaDataBrowser::updatePeptideHits(PeptideIdentification, int).<br>
        Updates the tree depending of the protein significance threshold.<br>
        Only ProteinHits with a score superior or equal to the current threshold will be displayed.
    */
    void updateTree_();

protected:
    /// Pointer to MetaDataBrowser
    MetaDataBrowser * pidv_caller_;
    /// The id of the item in the tree
    int tree_id_;

    ///@name Edit fields and buttons
    //@{
    QLineEdit * identifier_;
    QLineEdit * score_type_;
    QComboBox * higher_better_;
    QLineEdit * identification_threshold_;
    //@}

    /// Threshold for filtering by score
    QLineEdit * filter_threshold_;
  };
}
