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
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <algorithm>

namespace OpenMS
{

  TransformationDescription::TransformationDescription() :
    name_(), param_(), pairs_(), trafo_(0)
  {
  }

  TransformationDescription::~TransformationDescription()
  {
    delete trafo_;
  }

  TransformationDescription::TransformationDescription(const TransformationDescription& rhs) :
    name_(rhs.name_), param_(rhs.param_), pairs_(rhs.pairs_), trafo_(0)
  {
  }

  TransformationDescription& TransformationDescription::operator=(const TransformationDescription& rhs)
  {
    if (this == &rhs) return *this;

    name_ = rhs.name_;
    param_ = rhs.param_;
    pairs_ = rhs.pairs_;
		delete trafo_;
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

	const String& TransformationDescription::getName() const
	{
		return name_;
	}

  void TransformationDescription::setName(const String& name)
	{
		delete trafo_;
		trafo_ = 0;
		name_ = name;
	}
		
	const Param& TransformationDescription::getParameters() const
	{
		return param_;
	}
		
	void TransformationDescription::setParameters(const Param& param)
	{
		delete trafo_;
		trafo_ = 0;
		param_ = param;
	}
		
	const TransformationDescription::PairVector& TransformationDescription::getPairs() const
	{
		return pairs_;
	}

  TransformationDescription::PairVector& TransformationDescription::getPairs()
	{
		return pairs_;
	}
		
	void TransformationDescription::setPairs(const TransformationDescription::PairVector& pairs)
	{
		pairs_ = pairs;
	}

  const DataValue& TransformationDescription::getParam(const String& name) const
  {
    return param_.getValue(name);
  }

  void TransformationDescription::setParam(const String& name, DoubleReal value)
  {
    delete trafo_;
    trafo_ = 0;
    param_.setValue(name,value);
  }
		
	void TransformationDescription::setParam(const String& name, Int value)
	{
		delete trafo_;
		trafo_ = 0;
		param_.setValue(name,value);
	}
			
  void TransformationDescription::setParam(const String& name, const String& value)
  {
    delete trafo_;
    trafo_ = 0;
    param_.setValue(name,value);
  }

	void TransformationDescription::apply(DoubleReal& value) const
	{
		// initialize transformation (if unset).
		if (!trafo_) init_();
		// apply transformation
		trafo_->operator()(value);
	}

	void TransformationDescription::getInverse(TransformationDescription& result)
	{
		// initialize transformation (if unset).
		if (!trafo_) init_();
		// invert
		trafo_->getInverse(result);
	}

	DoubleReal TransformationDescription::getMaxRTErrorEstimate() const
	{
		// initialize transformation (if unset).
		if (!trafo_) init_();
		// invert
		return trafo_->getMaxRTErrorEstimate();
	}


//// internal structs

  struct TransformationDescription::None_ : TransformationDescription::Trafo_
  {
    None_(const TransformationDescription& rhs) :
      Trafo_(rhs)
    {
      return;
    }

    virtual void operator()(DoubleReal&) const
    {
      return;
    }

		virtual void getInverse(TransformationDescription& result)
		{
			if (this == result.trafo_) // self-assignment
			{
				result.pairs_.clear();
				result.param_.clear();
			}
			else
			{
				result.clear();
			}
			result.name_ = "none";
		}
		
  };

  struct TransformationDescription::Linear_ : TransformationDescription::Trafo_
  {
    Linear_(const TransformationDescription& rhs) :
      Trafo_(rhs)
    {
      if ( !rhs.param_.exists("slope") )
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "parameter 'slope' for 'linear' transformation not given");
      }
      slope_ = rhs.param_.getValue("slope");

      if ( !rhs.param_.exists("intercept") )
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "parameter 'intercept' for 'linear' transformation not given");
      }
      intercept_ = rhs.param_.getValue("intercept");
      return;
    }

    virtual void operator()(DoubleReal& value) const
    {
      value *= slope_;
      value += intercept_;
      return;
    }

		virtual void getInverse(TransformationDescription& result)
		{
			if (slope_ == 0) throw Exception::DivisionByZero(__FILE__, __LINE__, 
																											__PRETTY_FUNCTION__);

			if (this == result.trafo_) // self-assignment
			{
				result.pairs_.clear();
				result.param_.clear();
				intercept_ = -intercept_ / slope_;
				slope_ = 1.0 / slope_;
				result.param_.setValue("slope", slope_);
				result.param_.setValue("intercept", intercept_);
			}
			else
			{
				result.clear();
				result.param_.setValue("intercept", -intercept_ / slope_);
				result.param_.setValue("slope", 1.0 / slope_);
			}
			result.name_ = "linear";
		}

    protected:
      DoubleReal slope_;
      DoubleReal intercept_;
  };

  struct TransformationDescription::InterpolatedLinear_ : TransformationDescription::Trafo_
  {
    InterpolatedLinear_(const TransformationDescription& rhs) :
      Trafo_(rhs)
    {
      if ( rhs.pairs_.size() < 2 )
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "less than two pairs for 'interpolated_linear' transformation given");
      }
      pairs_ = rhs.pairs_;
      std::sort(pairs_.begin(), pairs_.end());
      return;
    }

    virtual void operator()(DoubleReal& value) const
    {
      if ( value <= pairs_.front().first )
      {
        DoubleReal slope = (pairs_.back().second - pairs_.front().second) / (pairs_.back().first - pairs_.front().first);
        value = pairs_.front().second + (value - pairs_.front().first) * slope;
      }
      else if ( value >= pairs_.back().first )
      {
        DoubleReal slope = (pairs_.back().second - pairs_.front().second) / (pairs_.back().first - pairs_.front().first);
        value = pairs_.back().second + (value - pairs_.back().first) * slope;
      }
      else
      {
        PairVector::const_iterator right = std::lower_bound(pairs_.begin(), pairs_.end(), PairVector::value_type(value, -std::numeric_limits<
            PairVector::value_type::second_type>::max()));
        PairVector::const_iterator left = right;
        --left;
        DoubleReal slope = (right->second - left->second) / (right->first - left->first);
        value = left->second + (value - left->first) * slope;
      }
      return;
    }

		virtual void getInverse(TransformationDescription& result)
		{
			if (this == result.trafo_) // self-assignment
			{
				result.param_.clear();
				for (PairVector::iterator it = pairs_.begin(); it != pairs_.end(); ++it)
				{
					*it = std::make_pair(it->second, it->first); // swap the components
				}
				std::sort(pairs_.begin(), pairs_.end());
				result.pairs_ = pairs_;
			}
			else
			{
				result.clear();
				for (PairVector::iterator it = pairs_.begin(); it != pairs_.end(); ++it)
				{
					result.pairs_.push_back(std::make_pair(it->second, it->first));
				}
			}
			result.name_ = "interpolated_linear";
		}

    protected:
      PairVector pairs_;
  };

  struct TransformationDescription::BSpline_ : TransformationDescription::Trafo_
  {
    BSpline_(const TransformationDescription& rhs) :
      Trafo_(rhs),
      tmp_lin_(rhs)
    {
      if (rhs.pairs_.size() < 4) // TODO: check number
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "less than four pairs for 'b_spline' transformation given (need at least four different pre-images)");
      }

      // TODO: allow arbitrary (non-uniform) breakpoints
      if (!rhs.param_.exists("num_breakpoints"))
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "parameter 'num_breakpoints' for 'b_spline' transformation missing");
      }

      size_ = rhs.pairs_.size();
      x_ = gsl_vector_alloc(size_);
      y_ = gsl_vector_alloc(size_);
      w_ = gsl_vector_alloc(size_);
      for (size_t i = 0; i < size_; ++i)
      {
        gsl_vector_set(x_, i, rhs.pairs_[i].first);
        gsl_vector_set(y_, i, rhs.pairs_[i].second);
        gsl_vector_set(w_, i, 1.0); // TODO: non-uniform weights
      }
      xmin_ = gsl_vector_min(x_);
      xmax_ = gsl_vector_max(x_);
		
      // set up cubic (k = 4) spline workspace:
      Size num_breakpoints = rhs.param_.getValue("num_breakpoints");
			if (num_breakpoints < 2) num_breakpoints = 2;
			else if (num_breakpoints > size_ - 2) num_breakpoints = size_ - 2;
			workspace_ = gsl_bspline_alloc(4, num_breakpoints);
			gsl_bspline_knots_uniform(xmin_, xmax_, workspace_);
      ncoeffs_ = gsl_bspline_ncoeffs(workspace_);
      xmin_ = gsl_vector_min(workspace_->knots);
      xmax_ = gsl_vector_max(workspace_->knots);
      computeFit_();
    }

    ~BSpline_()
    {
      gsl_bspline_free(workspace_);
      gsl_vector_free(bsplines_);
      gsl_vector_free(coeffs_);
      gsl_matrix_free(cov_);
      gsl_vector_free(x_);
      gsl_vector_free(y_);
      gsl_vector_free(w_);
    }

    void computeFit_()
    {
      // construct the fit matrix:
      gsl_matrix *fit_matrix = gsl_matrix_alloc(size_, ncoeffs_);
      bsplines_ = gsl_vector_alloc(ncoeffs_);
      for (size_t i = 0; i < size_; ++i)
      {
        double xi = gsl_vector_get(x_, i);
        gsl_bspline_eval(xi, bsplines_, workspace_);
        for (size_t j = 0; j < ncoeffs_; ++j)
        {
          double bspline = gsl_vector_get(bsplines_, j);
          gsl_matrix_set(fit_matrix, i, j, bspline);
        }
      }
      // do the fit:
      gsl_multifit_linear_workspace *multifit = gsl_multifit_linear_alloc(
				size_, ncoeffs_);
      coeffs_ = gsl_vector_alloc(ncoeffs_);
      cov_ = gsl_matrix_alloc(ncoeffs_, ncoeffs_);
      double chisq;
      gsl_multifit_wlinear(fit_matrix, w_, y_, coeffs_, cov_, &chisq, multifit);
      // clean-up:
      gsl_matrix_free(fit_matrix);
      gsl_multifit_linear_free(multifit);
      // for linear extrapolation (natural spline):
      computeLinear_(xmin_, slope_min_, offset_min_, sd_err_left_);
      computeLinear_(xmax_, slope_max_, offset_max_, sd_err_right_);
    }

    void computeLinear_(const double pos, double& slope, double& offset, double& sd_err)
    {
      gsl_bspline_deriv_workspace *deriv_workspace = gsl_bspline_deriv_alloc(4);
      gsl_matrix *deriv = gsl_matrix_alloc(ncoeffs_, 2);
      gsl_bspline_deriv_eval(pos, 1, deriv, workspace_, deriv_workspace);
      gsl_bspline_deriv_free(deriv_workspace);
      double results[2];
      for (size_t j = 0; j < 2; ++j)
      {
        for (size_t i = 0; i < ncoeffs_; ++i)
        {
          gsl_vector_set(bsplines_, i, gsl_matrix_get(deriv, i, j));
        }
        gsl_multifit_linear_est(bsplines_, coeffs_, cov_, &results[j], &sd_err);
        if (sd_err >= 1) // todo: find a good way to estimate the error!
        {
          LOG_ERROR << "BSpline::extrapolation() is unreliable. Consider reducing 'num_breakpoints' or use another model." << std::endl; 
        }
      }
      gsl_matrix_free(deriv);
      offset = results[0];
      slope = results[1];
    }

    virtual void operator()(DoubleReal& value) const
    {
      DoubleReal value_cpy = value;
      if (value < xmin_)
      { 
        value = offset_min_ - slope_min_ * (xmin_ - value);
        if (sd_err_left_ >= 1) bSplineCatch_(value, value_cpy);
      } 
      else if (value > xmax_) 
      { 
        value = offset_max_ + slope_max_ * (value - xmax_); 
        if (sd_err_right_ >= 1) bSplineCatch_(value, value_cpy);
      } 
      else // interpolation
      {
        double yerr;
        DoubleReal value_cpy = value;
        gsl_bspline_eval(value, bsplines_, workspace_);
        gsl_multifit_linear_est(bsplines_, coeffs_, cov_, &value, &yerr);
        if (yerr > 1) // the error gets too big!
        {
          bSplineCatch_(value, value_cpy);
        }
      }
    }

    /// alternative calculation of Trafo when BSpline gives bad SD
    /// , also updating 'max_rt_diff_'
    virtual void bSplineCatch_(const DoubleReal& value_bspline, DoubleReal& original_value) const
    {
      tmp_lin_(original_value);   // ... for linear interpolation
      max_rt_diff_ = std::max(max_rt_diff_, abs(original_value-value_bspline)); // get a feeling for the maximal error
    }

		virtual void getInverse(TransformationDescription& result)
		{
			if (this == result.trafo_) // self-assignment
			{
				result.param_.clear();
				result.pairs_.clear();
			}
			else
			{
				result.clear();
			}
			for (size_t i = 0; i < size_; ++i)
			{
				result.pairs_.push_back(std::make_pair(gsl_vector_get(y_, i),
																							 gsl_vector_get(x_, i)));
			}
			if (this == result.trafo_) // self-assignment
			{
				// swap x and y (weights stay the same):
				gsl_vector* temp;
				temp = x_;
				x_ = y_;
				y_ = temp;
				// recompute the spline:
				xmin_ = gsl_vector_min(x_);
				xmax_ = gsl_vector_max(x_);
				// workspace stays the same
				gsl_bspline_knots_uniform(xmin_, xmax_, workspace_);
				ncoeffs_ = gsl_bspline_ncoeffs(workspace_);
				xmin_ = gsl_vector_min(workspace_->knots);
				xmax_ = gsl_vector_max(workspace_->knots);
				computeFit_();
			}
			// ncoeffs = nbreak + k - 2
			result.param_.setValue("num_breakpoints", ncoeffs_ - 2);
			result.name_ = "b_spline";
		}

  protected:
    gsl_vector *x_, *y_, *w_, *bsplines_, *coeffs_;
    gsl_matrix *cov_;
    gsl_bspline_workspace *workspace_;
    size_t size_, ncoeffs_;
    double xmin_, xmax_; // first/last breakpoint
    double slope_min_, slope_max_, offset_min_, offset_max_;
    InterpolatedLinear_ tmp_lin_; // for extrapolation
    double sd_err_left_, sd_err_right_;
  };

  void TransformationDescription::init_() const
  {
    // workaround: init_() is const, but in fact it changes "hidden" state.
    Trafo_ * & trafo = const_cast<Trafo_*&> (trafo_);

    if (trafo) return; // object already present, delete it beforehand if you want a new one

    if (name_ == "none")
    {
      trafo = new None_(*this);
    }
    else if (name_ == "linear")
    {
      trafo = new Linear_(*this);
    }
    else if (name_ == "interpolated_linear")
    {
      trafo = new InterpolatedLinear_(*this);
    }
    else if (name_ == "b_spline")
    {
      trafo = new BSpline_(*this);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, (String("unknown transformation name '") + name_ + "'").c_str());
    }
  }

  std::ostream& operator<<(std::ostream& os, const TransformationDescription& td)
  {
    return os << " -- TransformationDescription  BEGIN --\n"
      "name: " << td.getName() << "\n"
      "parameters:\n" << td.getParameters() << " -- TransformationDescription END --" << std::endl;
  }

} // end of namespace OpenMS

