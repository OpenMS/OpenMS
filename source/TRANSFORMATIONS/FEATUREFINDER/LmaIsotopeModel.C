// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//									 OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//	Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//	This library is free software; you can redistribute it and/or
//	modify it under the terms of the GNU Lesser General Public
//	License as published by the Free Software Foundation; either
//	version 2.1 of the License, or (at your option) any later version.
//
//	This library is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
//	Lesser General Public License for more details.
//
//	You should have received a copy of the GNU Lesser General Public
//	License along with this library; if not, write to the Free Software
//	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA	02111-1307	USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaIsotopeModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <numeric>

namespace OpenMS
{
	LmaIsotopeModel::LmaIsotopeModel()
		: InterpolationModel(),
			charge_(0),
			monoisotopic_mz_(0.0)
	{
		setName(getProductName());

		defaults_.setValue("averagines:C",0.04443989f,"Number of C atoms per Dalton of mass.", StringList::create("advanced"));
		defaults_.setValue("averagines:H",0.06981572f,"Number of H atoms per Dalton of mass.", StringList::create("advanced"));
		defaults_.setValue("averagines:N",0.01221773f,"Number of N atoms per Dalton of mass.", StringList::create("advanced"));
		defaults_.setValue("averagines:O",0.01329399f,"Number of O atoms per Dalton of mass.", StringList::create("advanced"));
		defaults_.setValue("averagines:S",0.00037525f,"Number of S atoms per Dalton of mass.", StringList::create("advanced"));
		defaults_.setValue("isotope:trim_right_cutoff",0.001,"Cutoff in averagine distribution, trailing isotopes below this relative intensity are not considered.", StringList::create("advanced"));
		defaults_.setValue("isotope:maximum",100,"Maximum isotopic rank to be considered.", StringList::create("advanced"));
		defaults_.setValue("isotope:distance",1.000495,"Distance between consecutive isotopic peaks.", StringList::create("advanced"));
		defaults_.setValue("isotope:stdev",0.1,"Standard deviation of gaussian applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.", StringList::create("advanced"));
		defaults_.setValue("charge",1,"Charge state of the model.", StringList::create("advanced"));
		defaults_.setValue("statistics:mean",0.0,"Centroid m/z (as opposed to monoisotopic m/z).", StringList::create("advanced"));
		defaults_.setValue("bounding_box:min",0.0,"Lower end of bounding box enclosing the data used to fit the model.", StringList::create("advanced"));
		defaults_.setValue("bounding_box:max",1.0,"Upper end of bounding box enclosing the data used to fit the model.", StringList::create("advanced"));
		defaults_.setValue("total_intensity",100.0,"Total intensity under isotope curve.", StringList::create("advanced"));
		defaults_.setValue("monoisotopic_mz",0.0,"Position (mz) of the monoisotopic peak.", StringList::create("advanced"));

		defaultsToParam_();
	}

	LmaIsotopeModel::LmaIsotopeModel(const LmaIsotopeModel& source)
		: InterpolationModel(source)
	{
		setParameters( source.getParameters() );
		updateMembers_();
	}

	LmaIsotopeModel::~LmaIsotopeModel()
	{
	}

	LmaIsotopeModel& LmaIsotopeModel::operator = (const LmaIsotopeModel& source)
	{
		if (&source == this) return *this;

		InterpolationModel::operator = (source);
		setParameters( source.getParameters() );
		updateMembers_();

		return *this;
	}

	void LmaIsotopeModel::setSamples()
	{
		LinearInterpolation::container_type& data = interpolation_.getData();
		data.clear();
		if (max_==min_) return;
		data.reserve( UInt ( (max_-min_) / interpolation_step_ + 1 ) );
		CoordinateType pos = min_;

		// compute relative abundance of i-th isotopic peak
		LinearInterpolation::container_type isotopes_exact;
		CoordinateType mass = mean_ * charge_;

		Int C_num = Int( 0.5 + mass * averagine_[C]);
		Int N_num = Int( 0.5 + mass * averagine_[N]);
		Int O_num = Int( 0.5 + mass * averagine_[O]);
		Int H_num = Int( 0.5 + mass * averagine_[H]);
		Int S_num = Int( 0.5 + mass * averagine_[S]);

		String form("");
		if (C_num) form.append("C").append(String(C_num));
		if (H_num) form.append("H").append(String(H_num));
		if (N_num) form.append("N").append(String(N_num));
		if (O_num) form.append("O").append(String(O_num));
		if (S_num) form.append("S").append(String(S_num));

		// compute relative abundance of i-th isotopic peak of a peptide
		EmpiricalFormula formula(form);
		typedef IsotopeDistribution::iterator IsoIter;
		IsotopeDistribution isotope_distribution = formula.getIsotopeDistribution(max_isotope_);
		isotope_distribution.trimRight(trim_right_cutoff_);
		isotope_distribution.renormalize();
		for (IsoIter iter = isotope_distribution.begin(); iter != isotope_distribution.end(); ++iter)
		{
			isotopes_exact.push_back(iter->second);
		}

		CoordinateType term1 = 0;
		CoordinateType termSum = 0;
		for (Size step = 0; pos< max_; ++step)
		{
			// next pos
			pos = min_ + step * interpolation_step_;

			term1 = total_intensity_/(sqrt(2*Constants::PI)*isotope_stdev_);
			termSum = 0;
			for (Size i=0; i < isotopes_exact.size(); ++i)
			{
				termSum += isotopes_exact[i]*exp(-pow(pos-monoisotopic_mz_-i*isotope_distance_,2)/(2*isotope_stdev_*isotope_stdev_));
			}

			data.push_back(term1*termSum);
		}

		interpolation_.setScale	 ( interpolation_step_ );
		interpolation_.setOffset ( min_ );
	}

	void LmaIsotopeModel::setOffset(CoordinateType offset)
	{
		DoubleReal diff = offset - getInterpolation().getOffset();
		min_ += diff;
		max_ += diff;
		mean_ += diff;
		monoisotopic_mz_ += diff;

		InterpolationModel::setOffset(offset);

		param_.setValue("bounding_box:min", min_);
		param_.setValue("bounding_box:max", max_);
		param_.setValue("statistics:mean", mean_);
	}

	LmaIsotopeModel::CoordinateType LmaIsotopeModel::getOffset()
	{
		return getInterpolation().getOffset();
	}

	UInt LmaIsotopeModel::getCharge()
	{
		return charge_;
	}

	LmaIsotopeModel::CoordinateType LmaIsotopeModel::getCenter() const
	{
		return monoisotopic_mz_;
	}

	void LmaIsotopeModel::updateMembers_()
	{
		InterpolationModel::updateMembers_();

		monoisotopic_mz_ = param_.getValue("monoisotopic_mz");
		charge_ = param_.getValue("charge");
		isotope_stdev_ = param_.getValue("isotope:stdev");
		mean_ = param_.getValue("statistics:mean");
		max_isotope_ = param_.getValue("isotope:maximum");
		trim_right_cutoff_ = param_.getValue("isotope:trim_right_cutoff");
		isotope_distance_ = param_.getValue("isotope:distance");

		min_ = param_.getValue("bounding_box:min");
		max_ = param_.getValue("bounding_box:max");
		total_intensity_ = param_.getValue("total_intensity");

		averagine_[C] = param_.getValue("averagines:C");
		averagine_[H] = param_.getValue("averagines:H");
		averagine_[N] = param_.getValue("averagines:N");
		averagine_[O] = param_.getValue("averagines:O");
		averagine_[S] = param_.getValue("averagines:S");

		setSamples();
	}
}
