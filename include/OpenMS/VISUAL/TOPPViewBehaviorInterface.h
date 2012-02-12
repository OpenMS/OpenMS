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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPVIEWBEHAVIORINTERFACE_H
#define OPENMS_VISUAL_TOPPVIEWBEHAVIORINTERFACE_H

#include <QtCore/QObject>

namespace OpenMS
{
  /** @brief Interface class to model different behaviors of TOPPView.

   @ingroup TOPPView
   */
  class TOPPViewBehaviorInterface :
      public QObject
  {
    Q_OBJECT

  public:
    /// Destructor
    virtual ~TOPPViewBehaviorInterface() {};

  public slots:
    /// Behavior for showSpectraumAs1D
    virtual void showSpectrumAs1D(int index) = 0;

    /// Behavior for activate1DSpectrum
    virtual void activate1DSpectrum(int index) = 0;

    /// Behavior for deactivate1DSpectrum
    virtual void deactivate1DSpectrum(int index) = 0;

    /// Slot for behavior activation
    virtual void activateBehavior() = 0;

    /// Slot for behavior deactivation
    virtual void deactivateBehavior() = 0;
  };
}

#endif // OPENMS_VISUAL_TOPPVIEWBEHAVIORINTERFACE_H
