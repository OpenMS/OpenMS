// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_VISUALIZER_BASEVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_BASEVISUALIZER_H

#include <OpenMS/config.h>

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
    void load(ObjectType & o)
    {
      ptr_ = &o;
      temp_ = o;

      update_();
    }

protected:

    /// Pointer to the object that is currently edited
    ObjectType * ptr_;
    /// Copy of current object used to restore the original values
    ObjectType temp_;

protected:

    ///Updates the GUI from the temp_ variable.
    virtual void update_()
    {
    }

  };


}
#endif
