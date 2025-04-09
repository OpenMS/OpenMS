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
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

//STL
#include <vector>
#include <utility>

class QAbstractButton;
class QButtonGroup;

namespace OpenMS
{
  /**
      @brief MetaInfoVisualizer is a visualizer class for all classes that use one MetaInfo object as member.

      Meta information is an array of Type-Name-Value tuples. Classes that have a MetaInfo objects as a member can use this class to edit the MetaInfo object.
  */
  class OPENMS_GUI_DLLAPI MetaInfoVisualizer :
    public BaseVisualizerGUI,
    public BaseVisualizer<MetaInfoInterface>
  {
    Q_OBJECT

public:
    ///Constructor
    MetaInfoVisualizer(bool editable = false, QWidget * parent = nullptr);

    //Docu in base class
    void load(MetaInfoInterface & m);

public slots:

    //Docu in base class
    void store() override;

protected slots:

    /// Adds a new Type-Value pair to the MetaInfo Object.
    void add_();
    /// Clears out all fields.
    void clear_();
    /// Removes a selected Type-Value pair from the MetaInfo Object.
    void remove_(int);
    ///Undo the changes made in the GUI.
    void undo_();

protected:
    /// Loads all Type-Value pairs one after another.
    void loadData_(UInt index);

    /** @name Edit fields for new Type-Value pair.
      */
    //@{
    QLineEdit * newkey_;
    QLineEdit * newvalue_;
    QLineEdit * newdescription_;
    //@}

    ///@name Arrays of pointers to objects for temporary metaInfo data
    //@{
    std::vector<std::pair<UInt, QLineEdit *> > metainfoptr_;
    std::vector<std::pair<UInt, QLabel *> > metalabels_;
    std::vector<std::pair<UInt, QAbstractButton *> > metabuttons_;
    //@}

    ///@name Edit fields and buttons
    //@{
    QPushButton * addbutton_;
    QPushButton * clearbutton_;
    QButtonGroup * buttongroup_;
    //@}

    /// Counter to keep track of the actual row in the layout.
    int nextrow_;

    /// The layout to display the Type-Value pairs.
    QGridLayout * viewlayout_;

    /// Container for metainfo data.
    std::vector<UInt> keys_;
  };


}
