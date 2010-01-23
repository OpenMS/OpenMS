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
// $Maintainer: Stephan Aiche$
// $Authors: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_ISOTOPEMODELGENERAL_H
#define OPENMS_SIMULATION_ISOTOPEMODELGENERAL_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>

#include <numeric> // accumulate

namespace OpenMS
{

	/**
 		@brief A "general" isotope model which derives from class OpenMS::IsotopeModel.

		In contrast to the feature detection algorithm in OpenMS (from which we borrow this class),
		we don't need to calculate an average isotope distribution for a given mass but have an empirical formula
		and need its exact distribution.
 	*/
  class OPENMS_DLLAPI IsotopeModelGeneral
    : public IsotopeModel
  {
  public:

    /// Default constructor
    IsotopeModelGeneral();

    ///  copy constructor
    IsotopeModelGeneral(const IsotopeModelGeneral& source);

    /// destructor
    virtual ~IsotopeModelGeneral();

    /// assignment operator
    virtual IsotopeModelGeneral& operator =(const IsotopeModelGeneral& source);

		/// return charge state 
    //UInt getCharge();

    /// create new IsotopeModel object (needed by Factory)
    static BaseModel<1>* create()
    {
      return new IsotopeModelGeneral();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "IsotopeModelGeneral";
    }

    /// Initializes the model for a given formula. The formula should include any adducts (usually H+) that account for charge.
    void setSamples(EmpiricalFormula formula);

  protected:

    void updateMembers_();

  };

} // namespace OpenMS

#endif // OPENMS_SIMULATION_ISOTOPEMODELGENERAL_H
