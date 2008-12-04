// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

namespace OpenMS
{
	
	TransformationDescription::TransformationDescription()
		: name_(),
		  param_(),
		  pairs_(),
		  trafo_(0)
	{
	}
		
  TransformationDescription::~TransformationDescription()
  {
  	delete trafo_;
  }

	TransformationDescription::TransformationDescription(const TransformationDescription& rhs)
		: name_(rhs.name_),
			param_(rhs.param_),
		  pairs_(rhs.pairs_),
			trafo_(0)
	{
	}

  TransformationDescription& TransformationDescription::operator = (const TransformationDescription& rhs)
  {
    if (this==&rhs) return *this;
    
		name_ = rhs.name_;
		param_ = rhs.param_;
		pairs_ = rhs.pairs_;
		trafo_ = 0;
		
    return *this;
  }

	void TransformationDescription::clear()
	{
		name_ = "";
		param_.clear();
		pairs_.clear();
		delete trafo_;
		trafo_ = 0;
	}

	struct TransformationDescription::None_ : TransformationDescription::Trafo_
	{
		None_(const TransformationDescription& rhs)
			: Trafo_(rhs)
		{
			return;
		}

		virtual void operator ()(DoubleReal& ) const
		{
			return;
		}
	};
	
	struct TransformationDescription::Linear_ : TransformationDescription::Trafo_
	{
		Linear_(const TransformationDescription& rhs)
			: Trafo_(rhs)
		{
			if (!rhs.param_.exists("slope"))
			{
				throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"parameter 'slope' for 'linear' transformation not given");
			}
			slope_ = rhs.param_.getValue("slope");

			if (!rhs.param_.exists("intercept"))
			{
				throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"parameter 'intercept' for 'linear' transformation not given");
			}
			intercept_ = rhs.param_.getValue("intercept");
			return;
		}

		virtual void operator ()(DoubleReal& value) const
		{
			value *= slope_;
			value += intercept_;
			return;
		}

	 protected:
		DoubleReal slope_;
		DoubleReal intercept_;
	};
		
	struct TransformationDescription::InterpolatedLinear_ : TransformationDescription::Trafo_
	{
		InterpolatedLinear_(const TransformationDescription& rhs)
			: Trafo_(rhs)
		{
			if ( rhs.pairs_.size() < 2 )
			{
				throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"less than two pairs for 'interpolated_linear' transformation given");
			}
			pairs_ = rhs.pairs_;
			std::sort(pairs_.begin(),pairs_.end());
			return;
		}
		
		virtual void operator() (DoubleReal& value) const
		{
			if ( value <= pairs_.front().first )
			{
				DoubleReal slope = ( pairs_.back().second - pairs_.front().second ) / ( pairs_.back().first - pairs_.front().first );
				value = pairs_.front().second + ( value - pairs_.front().first ) * slope;
			}
			else if ( value >= pairs_.back().first )
			{
				DoubleReal slope = ( pairs_.back().second - pairs_.front().second ) / ( pairs_.back().first - pairs_.front().first );
				value = pairs_.back().second + ( value - pairs_.back().first ) * slope;
			}
			else
			{
				PairVector::const_iterator right =
					std::lower_bound( pairs_.begin(),
														pairs_.end(),
														PairVector::value_type( value,
																										- std::numeric_limits<PairVector::value_type::second_type>::max()
																									)
													);
				PairVector::const_iterator left = right;
				--left;
				DoubleReal slope = ( right->second - left->second ) / ( right->first - left->first );
				value = left->second + ( value - left->first ) * slope;
			}
			return;
		}

	 protected:
		PairVector pairs_;
	};

	void TransformationDescription::init_() const
	{
		// workaround: init_() is const, but in fact it changes "hidden" state.
		Trafo_ * & trafo = const_cast<Trafo_*&>(trafo_);

		if ( trafo ) delete trafo;
		trafo = 0;
		if (name_=="none")
		{
			trafo = new None_(*this);
		}
		else if (name_=="linear")
		{
			trafo = new Linear_(*this);
		}
		else if (name_=="interpolated_linear")
		{
			trafo = new InterpolatedLinear_(*this);
		}
		else
		{
			throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,(String("unknown transformation name '") + name_ + "'").c_str());
		}
	}




	std::ostream& operator<<(std::ostream& os, TransformationDescription const & td)
	{
		return os <<
		" -- TransformationDescription  BEGIN --\n"
		"name: " << td.getName() << "\n"
		"parameters:\n" << td.getParameters() <<
		" -- TransformationDescription END --" <<
		std::endl;
	}

} // end of namespace OpenMS

