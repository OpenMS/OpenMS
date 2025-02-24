// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>
#include <OpenMS/METADATA/Gradient.h>

//STL
#include <vector>

class QIntValidator;

namespace OpenMS
{
  /**
      @brief GradientVisualizer is a visualizer class for objects of type gradient.

      Each HPLC objects contains a gradient object. A gradient objects contains a list of eluents, timepoints and percentage values. Values can be added to the list, or the whole list can be deleted.
  */
  class OPENMS_GUI_DLLAPI GradientVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<Gradient>
  {
    Q_OBJECT

public:

    ///Constructor
    GradientVisualizer(bool editable = false, QWidget * parent = nullptr);

    //Docu in base class
    void load(Gradient & g);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    /// Add new timepoint to the list
    void addTimepoint_();
    /// Add new eluent to the list
    void addEluent_();
    ///Delete all data from gradient
    void deleteData_();
    ///Undo the changes made in the GUI.
    void undo_();

protected:
    /// Loads a list of eluent, timepoint and percentage triplets.
    void loadData_();
    /// Remove all data from layout
    void removeData_();


    /** @name Edit fields for new eluent-timepoint-percentage-triplets.
      */
    //@{
    QLineEdit * new_eluent_;
    QLineEdit * new_timepoint_;
    //@}

    /** @name Arrays of string values containing eluent, timepoint and percentage values.
    */
    //@{
    std::vector<String> eluents_;
    std::vector<Int> timepoints_;
    //@}

    /** @name Some buttons.
    */
    //@{
    QPushButton * add_eluent_button_;
    QPushButton * add_timepoint_button_;
    QPushButton * removebutton_;
    //@}

    /// Array of temporary pointers to gradient edit fields
    std::vector<QLineEdit *> gradientdata_;

    /// Array of temporary pointers to gradient labels
    std::vector<QLabel *> gradientlabel_;

    /// Pointer to fields with actual data
    QLineEdit * percentage_;

    /// A validator to check the input for the new timepoint.
    QIntValidator * timepoint_vali_;

    /// Counter to keep track of the actual row in the layout.
    int nextrow_;

    /// The layout to display the eluents, timepoints and percentages.
    QGridLayout * viewlayout_;

    //Docu in base class
    void update_() override;
  };


}
