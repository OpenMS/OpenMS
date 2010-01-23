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

  TransformationDescription&
  TransformationDescription::operator =(const TransformationDescription& rhs)
  {
    if ( this == &rhs )
      return *this;

    name_ = rhs.name_;
    param_ = rhs.param_;
    pairs_ = rhs.pairs_;
    trafo_ = 0;

    return *this;
  }

  void
  TransformationDescription::clear()
  {
    name_ = "";
    param_.clear();
    pairs_.clear();
    delete trafo_;
    trafo_ = 0;
  }

  struct TransformationDescription::None_ : TransformationDescription::Trafo_
  {
    None_(const TransformationDescription& rhs) :
      Trafo_(rhs)
    {
      return;
    }

    virtual void
    operator ()(DoubleReal&) const
    {
      return;
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

    virtual void
    operator ()(DoubleReal& value) const
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

    virtual void
    operator()(DoubleReal& value) const
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

  protected:
    PairVector pairs_;
  };

  struct TransformationDescription::BSpline_ : TransformationDescription::Trafo_
  {
    BSpline_(const TransformationDescription& rhs) :
      Trafo_(rhs)
    {
      if ( rhs.pairs_.size() < 4 ) // TODO: check number
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "less than four pairs for 'b_spline' transformation given (need at least four different pre-images)");
      }

      // TODO: allow arbitrary (non-uniform) breakpoints
      if ( !rhs.param_.exists("num_breakpoints") )
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "parameter 'num_breakpoints' for 'b_spline' transformation not given");
      }

      size = rhs.pairs_.size();
      x = gsl_vector_alloc(size);
      y = gsl_vector_alloc(size);
      w = gsl_vector_alloc(size);
      for ( size_t i = 0; i < size; ++i )
      {
        gsl_vector_set(x, i, rhs.pairs_[i].first);
        gsl_vector_set(y, i, rhs.pairs_[i].second);
        gsl_vector_set(w, i, 1.0); // TODO: non-uniform weights
      }
      xmin = gsl_vector_min(x);
      xmax = gsl_vector_max(x);
		
      // set up cubic (k = 4) spline workspace:
      Size num_breakpoints = rhs.param_.getValue("num_breakpoints");
			if (num_breakpoints < 2) num_breakpoints = 2;
			else if (num_breakpoints > size - 2) num_breakpoints = size - 2;
			workspace = gsl_bspline_alloc(4, num_breakpoints);
			gsl_bspline_knots_uniform(xmin, xmax, workspace);
      ncoeffs = gsl_bspline_ncoeffs(workspace);
      xmin = gsl_vector_min(workspace->knots);
      xmax = gsl_vector_max(workspace->knots);
      compute_fit_();
    }

    ~BSpline_()
    {
      gsl_bspline_free(workspace);
      gsl_vector_free(bsplines);
      gsl_vector_free(coeffs);
      gsl_matrix_free(cov);
      gsl_vector_free(x);
      gsl_vector_free(y);
      gsl_vector_free(w);
    }

    void
    compute_fit_()
    {
      // construct the fit matrix:
      gsl_matrix *fit_matrix = gsl_matrix_alloc(size, ncoeffs);
      bsplines = gsl_vector_alloc(ncoeffs);
      for ( size_t i = 0; i < size; ++i )
      {
        double xi = gsl_vector_get(x, i);
        gsl_bspline_eval(xi, bsplines, workspace);
        for ( size_t j = 0; j < ncoeffs; ++j )
        {
          double bspline = gsl_vector_get(bsplines, j);
          gsl_matrix_set(fit_matrix, i, j, bspline);
        }
      }
      // do the fit:
      gsl_multifit_linear_workspace *multifit = gsl_multifit_linear_alloc(
				size, ncoeffs);
      coeffs = gsl_vector_alloc(ncoeffs);
      cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
      double chisq;
      gsl_multifit_wlinear(fit_matrix, w, y, coeffs, cov, &chisq, multifit);
      // clean-up:
      gsl_matrix_free(fit_matrix);
      gsl_multifit_linear_free(multifit);
      // for linear extrapolation (natural spline):
      compute_linear_(xmin, slope_min, offset_min);
      compute_linear_(xmax, slope_max, offset_max);
    }

    void
    compute_linear_(double pos, double& slope, double& offset)
    {
      gsl_bspline_deriv_workspace *deriv_workspace = gsl_bspline_deriv_alloc(4);
      gsl_matrix *deriv = gsl_matrix_alloc(ncoeffs, 2);
      gsl_bspline_deriv_eval(pos, 1, deriv, workspace, deriv_workspace);
      gsl_bspline_deriv_free(deriv_workspace);
      double results[2];
      for ( size_t j = 0; j < 2; ++j )
      {
        for ( size_t i = 0; i < ncoeffs; ++i )
        {
          gsl_vector_set(bsplines, i, gsl_matrix_get(deriv, i, j));
        }
        double yerr;
        gsl_multifit_linear_est(bsplines, coeffs, cov, &results[j], &yerr);
      }
      gsl_matrix_free(deriv);
      offset = results[0];
      slope = results[1];
    }

    virtual void
    operator()(DoubleReal& value) const
    {
      if ( value < xmin )
      {
        value = offset_min - slope_min * (xmin - value);
      }
      else if ( value > xmax )
      {
        value = offset_max + slope_max * (value - xmax);
      }
      else
      {
        double yerr;
        gsl_bspline_eval(value, bsplines, workspace);
        gsl_multifit_linear_est(bsplines, coeffs, cov, &value, &yerr);
      }
    }

  protected:
    gsl_vector *x, *y, *w, *bsplines, *coeffs;
    gsl_matrix *cov;
    gsl_bspline_workspace *workspace;
    size_t size, ncoeffs;
    double xmin, xmax; // first/last breakpoint
    double slope_min, slope_max, offset_min, offset_max;
  };

  void
  TransformationDescription::init_() const
  {
    // workaround: init_() is const, but in fact it changes "hidden" state.
    Trafo_ * & trafo = const_cast<Trafo_*&> (trafo_);

    if ( trafo )
      delete trafo;
    trafo = 0;
    if ( name_ == "none" )
    {
      trafo = new None_(*this);
    }
    else if ( name_ == "linear" )
    {
      trafo = new Linear_(*this);
    }
    else if ( name_ == "interpolated_linear" )
    {
      trafo = new InterpolatedLinear_(*this);
    }
    else if ( name_ == "b_spline" )
    {
      trafo = new BSpline_(*this);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, (String("unknown transformation name '") + name_ + "'").c_str());
    }
  }

  std::ostream&
  operator<<(std::ostream& os, TransformationDescription const & td)
  {
    return os << " -- TransformationDescription  BEGIN --\n"
      "name: " << td.getName() << "\n"
      "parameters:\n" << td.getParameters() << " -- TransformationDescription END --" << std::endl;
  }

} // end of namespace OpenMS

