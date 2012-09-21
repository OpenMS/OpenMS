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
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdio.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Info.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FeatureLCProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SHFeature.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMS.h>

namespace OpenMS
{

	using namespace std;

	double LCMS::MINIMAL_PEP_PROPHET_THERSHOLD = -3.0;
//	double LCMS::PEP_PROPHET_THERSHOLD;

////////////////////////////////////////////////
// constructor 
	LCMS::LCMS()
	{
		spectrum_id = -1;
		MASTER_ID = -1;
	}

////////////////////////////////////////////////
// constructor for the object LCMS:
	LCMS::LCMS(string IN_name)
	{
		spec_name = IN_name;
		spectrum_id = -1;
		MASTER_ID = -1;
	}

//////////////////////////////////////////////////
// copy constructor:
	LCMS::LCMS(const LCMS* tmp)
	{
		spec_name = tmp->spec_name;
		spectrum_id = tmp->spectrum_id;
		raw_spec_names = tmp->raw_spec_names;
		MASTER_ID = tmp->MASTER_ID;
		ALIGNMENT_ERROR = tmp->ALIGNMENT_ERROR;
		feature_list = tmp->feature_list;
	}

//////////////////////////////////////////////////
// copy constructor:
	LCMS::LCMS(const LCMS& tmp)
	{
		spec_name = tmp.spec_name;
		spectrum_id = tmp.spectrum_id;
		raw_spec_names = tmp.raw_spec_names;
		MASTER_ID = tmp.MASTER_ID;
		ALIGNMENT_ERROR = tmp.ALIGNMENT_ERROR;
		feature_list = tmp.feature_list;
	}

//////////////////////////////////////////////////
// copy constructor:
	LCMS& LCMS::operator=(const LCMS& tmp)
	{

		spec_name = tmp.spec_name;
		spectrum_id = tmp.spectrum_id;
		raw_spec_names = tmp.raw_spec_names;
		MASTER_ID = tmp.MASTER_ID;
		ALIGNMENT_ERROR = tmp.ALIGNMENT_ERROR;
		feature_list = tmp.feature_list;

		return *this;
	}

//////////////////////////////////////////////////
// class desctructor
	LCMS::~LCMS()
	{

		//  destroy the feature list:
		if (!feature_list.empty())
			feature_list.clear();

		if (!raw_spec_names.empty())
			raw_spec_names.clear();

		ALIGNMENT_ERROR.clear();

	}

//////////////////////////////////////////////////
// order the features according to parent mass:
	void LCMS::order_by_mass()
	{
		sort(feature_list.begin(), feature_list.end(), OPERATOR_MZ());
	}

//////////////////////////////////////////////////
// show the content of the spectra:
	void LCMS::show_info()
	{

		if (get_spec_name().size() > 0)
			printf("\t\t -- LC-MS name: %s ", get_spec_name().c_str());
		else
			printf("\t\t -- LC-MS ID: %d,", get_spectrum_ID());

		if (get_nb_raw_specs() != 0)
		{
			printf("[MASTER MAP ID=%d] ", get_MASTER_ID());
		}
		else
		{
			printf("[LC-MS ID=%d] ", get_spectrum_ID());
		}

		printf(" #features: %d, #MS/MS ids: %d (no Thresholding: %d)\n", get_nb_features(), get_nb_identified_features(),
				get_nb_identified_features(MINIMAL_PEP_PROPHET_THERSHOLD));

		// show the child LC/MS runs:
		map<int, string>::iterator C = get_raw_spec_name_start();
		while (C != get_raw_spec_name_end())
		{
			printf("\t\t\t - Child LC-MS: %s [ID=%d]\n", (*C).second.c_str(), (*C).first);
			C++;
		}
		vector<SHFeature>::iterator p = feature_list.begin();
		while (p != feature_list.end())
		{
			// if((*p).get_MS2_info()){
			//(*p).show_info();
			//}
			p++;
		}
	}

//////////////////////////////////////////////////
// search the list of feature for the one with input ID:
	SHFeature* LCMS::find_feature_by_ID(int ID)
	{

		vector<SHFeature>::iterator p = feature_list.begin();
		while (p != feature_list.end())
		{

			if (ID == (*p).get_feature_ID())
			{
				return &(*p);
			}

			p++;
		}

		return NULL;
	}

//////////////////////////////////////////////////
// count the number of common peaks of a given number of LC-MS:
// ONLY count the one which appear count-times
	int LCMS::get_nb_common_peaks(int count)
	{

		// search for the others:
		int common_count = 0;

		SHFeature* PEAK = NULL;
		vector<SHFeature>::iterator p = feature_list.begin();
		while (p != feature_list.end())
		{

			// get the peak at a charge state:
			PEAK = &(*p);

			if (PEAK != NULL)
			{
				if (PEAK->get_nb_common_match() == count)
				{
					common_count++;
				}
			}
			// next feature
			p++;
		}

		return common_count;
	}

//////////////////////////////////////////////////////////////////
// remove a feature from the LC/MS run:
	void LCMS::remove_feature(SHFeature* IN)
	{
		vector<SHFeature>::iterator P = find(feature_list.begin(), feature_list.end(), IN);
		if (P != feature_list.end())
		{
			P->show_info();
			P = feature_list.erase(P);
		}
	}

//////////////////////////////////////////////////////////////////
// remove a feature from the LC/MS run by ID:
	void LCMS::remove_feature_by_ID(SHFeature* IN)
	{
		remove_feature_by_ID(IN->get_feature_ID());
	}

//////////////////////////////////////////////////////////////////
// remove a feature from the LC/MS run by ID:
	void LCMS::remove_feature_by_ID(int ID)
	{
		vector<SHFeature>::iterator P = feature_list.begin();
		while (P != feature_list.end())
		{
			if (P->get_feature_ID() == ID)
			{
				feature_list.erase(P);
				break;
			}

			P++;

		}
	}

///////////////////////////////////////////////////////////
// get alignment error at specific TR:
	void LCMS::get_alignment_error(double IN, double* UP, double* DOWN)
	{

		if (!ALIGNMENT_ERROR.empty())
		{

			map<double, pair<double, double> >::iterator P = ALIGNMENT_ERROR.lower_bound(IN);

			// exact match:
			if ((*P).first == IN)
			{
				*UP = (*P).second.first;
				*DOWN = (*P).second.second;
			}

			// is bigger than any element:
			else if (P == ALIGNMENT_ERROR.end())
			{
				P--;
				*UP = (*P).second.first;
				*DOWN = (*P).second.second;
			}

			// is smaller than any element:
			else if (P == ALIGNMENT_ERROR.begin())
			{
				*UP = (*P).second.first;
				*DOWN = (*P).second.second;
			}

			else
			{
				// is within the range of the list but has not found an exact match
				// -> extrapolate between 2 data points:
				double up_E_u = (*P).second.first;
				double down_E_u = (*P).second.second;
				double TR_u = (*P).first;

				P--;

				double up_E_d = (*P).second.first;
				double down_E_d = (*P).second.second;
				double TR_d = (*P).first;

				double w_u = (IN - TR_d) / (TR_u - TR_d);
				double w_d = (TR_u - IN) / (TR_u - TR_d);

				*UP = w_u * up_E_u + w_d * up_E_d;
				*DOWN = w_u * down_E_u + w_d * down_E_d;
			}

		}
	}

///////////////////////////////////////////////////////////////
// check if this LC/MS ID is present in the raw LC/MS runs 
	bool LCMS::find_LC_MS_by_ID(int ID)
	{

		map<int, string>::iterator F = raw_spec_names.find(ID);
		if (F != raw_spec_names.end())
		{
			return true;
		}
		else
		{
			return false;
		}
	}

///////////////////////////////////////////////////////////////
// compare the LC/MS runs names 
	bool LCMS::check_LCMS_name(string IN)
	{

		if (spec_name.find(IN) != string::npos)
		{
			return true;
		}

		map<int, string>::iterator F = raw_spec_names.begin();
		while (F != raw_spec_names.end())
		{
			if ((*F).second.find(IN) != string::npos)
			{
				return true;
			}
			F++;
		}

		return false;
	}

//////////////////////////////////////////////////
// set the id of all features
	void LCMS::setFeatureLCMSID()
	{

		vector<SHFeature>::iterator p = feature_list.begin();
		while (p != feature_list.end())
		{
			(*p).set_spectrum_ID(get_spectrum_ID());
			p++;
		}
	}

}
