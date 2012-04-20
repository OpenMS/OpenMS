// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENTHELPER_H
#define OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENTHELPER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/CVTerm.h>
#include <OpenMS/METADATA/CVTermList.h>

namespace OpenMS
{
  /**
    @brief This class stores helper structures that are used in multiple
    classes of the TargetedExperiment (e.g. ReactionMonitoringTransition and
    IncludeExcludeTarget).
  */

  struct OPENMS_DLLAPI TargetedExperimentHelper
  {

    struct Configuration
      : public CVTermList
    {
      String contact_ref;
      String instrument_ref;
      std::vector<CVTermList> validations;

      Configuration& operator = (const Configuration& rhs)
      {
        if (this != &rhs)
        {
          CVTermList::operator = (rhs);
          contact_ref = rhs.contact_ref;
          instrument_ref = rhs.instrument_ref;
          validations = rhs.validations;
        }
        return *this;
      }
    };

  };
}

#endif // OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENTHELPER_H

