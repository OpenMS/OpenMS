// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_HPLC_H
#define OPENMS_METADATA_HPLC_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Gradient.h>

namespace OpenMS 
{
  /**
    @brief Representation of a HPLC experiment
    
    It contains the description of instrument, the settings and the gradient.
		
		@ingroup Metadata
  */
  class OPENMS_DLLAPI HPLC
  {
		public:
      /// Constructor
      HPLC();
      /// Copy constructor
      HPLC(const HPLC& source);
      /// Destructor
      ~HPLC();
      
      /// Assignment operator
      HPLC& operator = (const HPLC& source);
      
      /// Equality operator
      bool operator == (const HPLC& source) const;
      /// Equality operator
      bool operator != (const HPLC& source) const;

      /// returns a const reference to the instument name
      const String& getInstrument() const;
      /// sets the instument name
      void setInstrument(const String& instrument);

      /// returns a const reference to the column description 
      const String& getColumn() const;
      /// sets the column description 
      void setColumn(const String& column);

      /// returns the temperature (in degree C)
      Int getTemperature() const;
      /// sets the temperature (in degree C)
      void setTemperature(Int temperature);

      /// returns the pressure (in bar)
      UInt getPressure() const;
      /// sets the pressure (in bar)
      void setPressure(UInt pressure);

      /// returns the flux (in microliter/sec)
      UInt getFlux() const;
      /// sets the flux (in microliter/sec)
      void setFlux(UInt flux);

      /// returns the comments
      String getComment() const;
      /// sets the comments
      void setComment(String comment);

      /// returns a const reference to the used gradient 
      const Gradient& getGradient() const;
      /// returns a mutable reference to the used gradient
      Gradient& getGradient();
      /// sets the used gradient
      void setGradient(const Gradient& gradient);

    protected:
    	String instrument_;
      String column_;
      Int temperature_;
      Int pressure_;
      Int flux_;
      String comment_;
      Gradient gradient_;
  };
 
} // namespace OpenMS

#endif // OPENMS_METADATA_HPLC_H


