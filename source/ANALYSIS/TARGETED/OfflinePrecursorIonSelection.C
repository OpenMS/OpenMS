// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/KERNEL/ComparatorUtils.h>

namespace OpenMS
{

OfflinePrecursorIonSelection::OfflinePrecursorIonSelection() : DefaultParamHandler("OfflinePrecursorIonSelection")
{  
	defaults_.setValue("ms2_spectra_per_rt_bin",5,"Number of allowed MS/MS spectra in a retention time bin.");
	defaults_.setMinInt("ms2_spectra_per_rt_bin",1);

	defaults_.setValue("min_peak_distance",3.,"The minimal distance (in Da) of two peaks in one spectrum so that they can be selected.");
	defaults_.setMinFloat("min_peak_distance",0.);

	defaults_.setValue("selection_window",2.,"All peaks within a mass window (in Da) of a selected peak are also selected for fragmentation.");
  defaults_.setMinFloat("selection_window",0.);

	defaults_.setValue("exclude_overlapping_peaks","false","If true overlapping or nearby peaks (within min_peak_distance) are excluded for selection.");
	defaults_.setValidStrings("exclude_overlapping_peaks",StringList::create("true,false"));

  defaults_.setValue("Exclusion:use_dynamic_exclusion","false","If true dynamic exclusion is applied.");
  defaults_.setValidStrings("Exclusion:use_dynamic_exclusion",StringList::create("true,false"));

  defaults_.setValue("Exclusion:exclusion_time",100.,"The time (in seconds) a feature is excluded.");
  defaults_.setMinFloat("Exclusion:exclusion_time",0.);
	defaultsToParam_();
}

OfflinePrecursorIonSelection::~OfflinePrecursorIonSelection()
{

}

}
