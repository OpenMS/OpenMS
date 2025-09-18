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

namespace OpenMS
{
  /**
      @brief A base class for all visualizer classes

      This class provides members and functions that depend on the visualizer type.

      The GUI components are provided by the BaseVisualizerGUI class.
      The two classes cannot be merged, as templates and the Qt meta object compiler cannot be combined.

      Visualizers are mainly used by the MetaDataBrowser.
  */
  template <typename ObjectType>
  class OPENMS_GUI_DLLAPI BaseVisualizer
  {
public:

    /// Loads the object that is to be edited.
    void load(ObjectType& o)
    {
      ptr_ = &o;
      temp_ = o;

      update_();
    }

    virtual ~BaseVisualizer() {}

protected:

    /// Pointer to the object that is currently edited
    ObjectType* ptr_;
    /// Copy of current object used to restore the original values
    ObjectType temp_;

protected:

    ///Updates the GUI from the temp_ variable.
    virtual void update_()
    {
    }

  };


}
