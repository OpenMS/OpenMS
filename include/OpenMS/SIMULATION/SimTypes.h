// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Stephan Aiche Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_SIMTYPES_H
#define OPENMS_SIMULATION_SIMTYPES_H

#include <vector>
#include <utility>
#include <map>

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

namespace OpenMS 
{
  /// Coordinate type in mz and rt dimension
  typedef Peak2D::CoordinateType SimCoordinateType;
  
	/// Abundance of proteins/peptides
	typedef Peak2D::IntensityType SimIntensityType;
	
	/// Sequence -> Intensity container
	typedef std::map< String, SimIntensityType > SampleProteins;
	/// Peptides and Proteins are the same structurewise
	typedef SampleProteins SamplePeptides;

  /// A posttranslational modification
  struct PTM
  {
    /// (Simplified) name
    String name_;
    
    /// Formula
    EmpiricalFormula formula_;
    
    /// Relative abundance (in %)
    double abundance_;
    
    /// mass shift (positive = true, negative = false)
    bool shift_;
  };  
  
  // maps from aminoacid (e.g. "A") to a list of possible modifications
  typedef std::multimap<String, PTM> PTMTable;
}

#endif
