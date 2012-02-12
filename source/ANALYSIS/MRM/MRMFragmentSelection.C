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

#include <OpenMS/ANALYSIS/MRM/MRMFragmentSelection.h>

#include <algorithm>

using namespace std;


namespace OpenMS
{
	MRMFragmentSelection::MRMFragmentSelection()
		: DefaultParamHandler("MRMFragmentSelection")
	{
		defaults_.setValue("num_top_peaks", 4, "Number of most intense peak to pick");
		defaults_.setValue("min_pos_precursor_percentage", 80.0, "Minimal ion position the ion should have, relative to the precursor position");
		defaults_.setValue("min_mz", 400.0, "Minimal m/z value that is allowed for selection.");
		defaults_.setValue("max_mz", 1200.0, "Maximal m/z value that is allowed for selection.");
		defaults_.setValue("consider_names", "true", "Should names be considered when selecting ions?");
		defaults_.setValidStrings("consider_names", StringList::create("true,false"));
		defaults_.setValue("allow_loss_ions", "false", "Should loss ions allowed to be selected?");
		defaults_.setValidStrings("allow_loss_ions", StringList::create("true,false"));
		defaults_.setValue("allowed_ion_types", StringList::create("y"), "The one-character-typenames of the ion types allowed"); 
		defaults_.setValue("allowed_charges", StringList::create("1"), "List of allowed charge states for selection.");

		defaultsToParam_();
	}

  MRMFragmentSelection::MRMFragmentSelection(const MRMFragmentSelection& rhs)
		: DefaultParamHandler(rhs)
	{
	}

	MRMFragmentSelection::~MRMFragmentSelection()
	{
	}

	MRMFragmentSelection& MRMFragmentSelection::operator = (const MRMFragmentSelection& rhs)
	{
		if (&rhs != this)
		{
			DefaultParamHandler::operator = (rhs);
		}
		return *this;
	}

	void MRMFragmentSelection::selectFragments(vector<RichPeak1D>& selected_peaks, const RichPeakSpectrum& spec)
	{
    Size num_top_peaks = param_.getValue("num_top_peaks");
		bool consider_names(param_.getValue("consider_names").toBool());
    DoubleReal min_pos_precursor_percentage = (DoubleReal)param_.getValue("min_pos_precursor_percentage") / 100.0;
		DoubleReal min_mz = (DoubleReal)param_.getValue("min_mz");
		DoubleReal max_mz = (DoubleReal)param_.getValue("max_mz");
		if (spec.getPrecursors().empty())
		{
			cerr << "MRMFragmentSelection: No Precursor peaks defined! Bailing out..." << endl;
			return;
		}
		DoubleReal precursor_pos =  spec.getPrecursors().begin()->getMZ();
    RichPeakSpectrum spec_copy = spec;
    spec_copy.sortByIntensity(true);

    for (Size i = 0; i < spec_copy.size() && selected_peaks.size() < num_top_peaks; ++i)
    {
      String name = spec_copy[i].getMetaValue("IonName");
      //if (spec_copy[i].metaValueExists("MSPPeakInfo"))
      //{
      //  name = spec_copy[i].getMetaValue("MSPPeakInfo");
      //}
      
			if (spec_copy[i].getMZ() >= min_mz && spec_copy[i].getMZ() <= max_mz &&
					spec_copy[i].getMZ() > min_pos_precursor_percentage * precursor_pos &&
         (!consider_names || peakselectionIsAllowed_(spec_copy[i])))
      {
        selected_peaks.push_back(spec_copy[i]);
      }
    }

    return;
	}

	bool MRMFragmentSelection::peakselectionIsAllowed_(const RichPeak1D& peak)
	{
		StringList allowed_charges = param_.getValue("allowed_charges");

		String name;
		if (peak.metaValueExists("IonName"))
		{
			name = peak.getMetaValue("IonName");
		}
		
		if (name != "")
		{
			StringList allowed_types((StringList)param_.getValue("allowed_ion_types"));
			bool type_found(false);
			for (StringList::const_iterator it = allowed_types.begin(); it != allowed_types.end(); ++it)
			{
				if (name.hasSubstring(*it))
				{
					type_found = true;
				}
			}
			if (type_found)
			{
				bool allow_loss_ions(param_.getValue("allow_loss_ions").toBool());
				Size charges = count(name.begin(), name.end(), '+');
				bool charges_ok = allowed_charges.contains(String(charges));
				if (allow_loss_ions && charges_ok)
				{
					// TODO implement charges
					return true;
				}
				else
				{
					if (! (name.hasSubstring("-H") || name.hasSubstring("-C") || name.hasSubstring("-N")))
					{
						Size c = count(name.begin(), name.end(), '+');
						if (allowed_charges.contains(c))
						{
							return true;
						}
						else
						{
							return false;
						}
					}
					else
					{
						return false;
					}
				}
			}
			else
			{
				return false;
			}
		}
		else
		{
			return false;
		}
		return true;
	}

}


