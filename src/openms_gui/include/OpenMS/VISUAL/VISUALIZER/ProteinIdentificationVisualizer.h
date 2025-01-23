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
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

namespace OpenMS
{
  class MetaDataBrowser;

  /**
      @brief Class that displays all meta information for ProteinIdentification objects

      This class provides all functionality to view the meta information of an object of type Identification.
  */
  class OPENMS_GUI_DLLAPI ProteinIdentificationVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<ProteinIdentification>
  {
    Q_OBJECT

public:
    ///Constructor
    ProteinIdentificationVisualizer(bool editable = false, QWidget * parent = 0, MetaDataBrowser * caller = nullptr);

    /// Loads the meta data from the object to the viewer. Gets the id of the item in the tree as parameter.
    void load(ProteinIdentification & s, int tree_item_id);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    /**
        @brief Updates the tree by calling MetaDataBrowser::updateProteinHits()

        Calls MetaDataBrowser::updateProteinHits().<br>
        Updates the tree depending of the protein significance threshold.<br>
        Only ProteinHits with a score superior or equal to the current threshold will be displayed.
    */
    void updateTree_();

    ///Undo the changes made in the GUI.
    void undo_();

protected:

    /// Pointer to MetaDataBrowser
    MetaDataBrowser * pidv_caller_;
    /// The id of the item in the tree
    int tree_id_;

    ///@name Edit fields and buttons
    //@{
    QLineEdit * engine_;
    QLineEdit * engine_version_;
    QLineEdit * identification_date_;
    QLineEdit * identification_threshold_;
    QLineEdit * identifier_;
    QLineEdit * score_type_;
    QComboBox * higher_better_;

    QLineEdit * db_;
    QLineEdit * db_version_;
    QLineEdit * taxonomy_;
    QLineEdit * charges_;
    QLineEdit * missed_cleavages_;
    QLineEdit * peak_tolerance_;
    QLineEdit * precursor_tolerance_;
    QComboBox * mass_type_;
    QLineEdit * enzyme_;
    //@}

    /// Threshold for filtering by score
    QLineEdit * filter_threshold_;
  };
}
