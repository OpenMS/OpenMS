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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaIsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <numeric>
#include <boost/math/special_functions/fpclassify.hpp>

namespace OpenMS
{

		LmaIsotopeFitter1D::LmaIsotopeFitter1D()
    : LevMarqFitter1D()
    {
      setName(getProductName());

      defaults_.setValue("averagines:C",0.0443f,"Number of C atoms per Dalton of mass.", StringList::create("advanced"));
      defaults_.setValue("averagines:H",0.007f,"Number of H atoms per Dalton of mass.", StringList::create("advanced"));
      defaults_.setValue("averagines:N",0.0037f,"Number of N atoms per Dalton of mass.", StringList::create("advanced"));
      defaults_.setValue("averagines:O",0.022f,"Number of O atoms per Dalton of mass.", StringList::create("advanced"));
      defaults_.setValue("averagines:S",0.00037f,"Number of S atoms per Dalton of mass.", StringList::create("advanced"));
      defaults_.setValue("isotope:trim_right_cutoff",0.001,"Cutoff in averagine distribution, trailing isotopes below this relative intensity are not considered.", StringList::create("advanced"));
      defaults_.setValue("isotope:maximum",100,"Maximum isotopic rank to be considered.", StringList::create("advanced"));
      defaults_.setValue("isotope:distance",1.000495,"Distance between consecutive isotopic peaks.", StringList::create("advanced"));
      defaults_.setValue("isotope:stdev",0.1,"Standard deviation of gaussian applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.", StringList::create("advanced"));
      defaults_.setValue("charge",1,"Charge state of the model.", StringList::create("advanced"));
      defaults_.setValue("statistics:mean",0.0,"Centroid m/z (as opposed to monoisotopic m/z).", StringList::create("advanced"));
      defaults_.setValue("statistics:variance",1.0,"Variance of the model.", StringList::create("advanced"));
      defaults_.setValue("interpolation_step",0.1,"Sampling rate for the interpolation of the model function.", StringList::create("advanced"));
      defaults_.setValue("total_intensity",100.0,"Total intensity under the curve in mz dimension.", StringList::create("advanced"));
      defaults_.setValue("monoisotopic_mass",0.0,"Monoisotopic mz of the model.", StringList::create("advanced"));

      defaultsToParam_();
    }

    LmaIsotopeFitter1D::LmaIsotopeFitter1D(const LmaIsotopeFitter1D& source)
    : LevMarqFitter1D(source)
    {
      updateMembers_();
    }

    LmaIsotopeFitter1D::~LmaIsotopeFitter1D()
    {
    }

    LmaIsotopeFitter1D& LmaIsotopeFitter1D::operator = (const LmaIsotopeFitter1D& source)
    {
      if (&source == this) return *this;

      LevMarqFitter1D::operator = (source);
      updateMembers_();

      return *this;
    }

    Int LmaIsotopeFitter1D::residual_(const gsl_vector* x, void* params, gsl_vector* f)
    {
        Size n = static_cast<LmaIsotopeFitter1D::Data*> (params) ->n;
        RawDataArrayType set = static_cast<LmaIsotopeFitter1D::Data*> (params) ->set;
        ContainerType isotopes_exact = static_cast<LmaIsotopeFitter1D::Data*> (params) ->isotopes_exact;
        CoordinateType isotope_distance = static_cast<LmaIsotopeFitter1D::Data*> (params) ->isotope_distance;
        CoordinateType sigma = static_cast<LmaIsotopeFitter1D::Data*> (params) ->sigma;
        //CoordinateType stdev = static_cast<LmaIsotopeFitter1D::Data*> (params) ->isotopes_stdev;

        CoordinateType A = gsl_vector_get( x, 0 );
        CoordinateType stdev = static_cast<LmaIsotopeFitter1D::Data*> (params) ->isotopes_stdev; //gsl_vector_get( x, 1 );
        CoordinateType mono_mz = gsl_vector_get( x, 1 );

        CoordinateType Yi = 0.0;

        // iterate over all points of the signal
        for ( Size i = 0; i < n; ++i )
        {
            CoordinateType m = set[i].getPos();

            CoordinateType term1 = A/(sqrt(2*Constants::PI)*stdev);
            CoordinateType termSum = 0;
            for (Size j=0; j<isotopes_exact.size(); ++j)
            {
              termSum += isotopes_exact[j]*exp(-pow(m-mono_mz-j*isotope_distance,2)/(2*stdev*stdev)) ;
            }

            Yi = term1*termSum;

            gsl_vector_set( f, i, (Yi - set[i].getIntensity())/sigma );
        }

        return GSL_SUCCESS;
    }

    Int LmaIsotopeFitter1D::jacobian_(const gsl_vector* x, void* params, gsl_matrix* J)
    {
        Size n = static_cast<LmaIsotopeFitter1D::Data*> (params) ->n;
        RawDataArrayType set = static_cast<LmaIsotopeFitter1D::Data*> (params) ->set;
        ContainerType isotopes_exact = static_cast<LmaIsotopeFitter1D::Data*> (params) ->isotopes_exact;
        CoordinateType isotope_distance = static_cast<LmaIsotopeFitter1D::Data*> (params) ->isotope_distance;
        CoordinateType sigma = static_cast<LmaIsotopeFitter1D::Data*> (params) ->sigma;
        //CoordinateType stdev = static_cast<LmaIsotopeFitter1D::Data*> (params) ->isotopes_stdev;

        CoordinateType A = gsl_vector_get( x, 0 );
        CoordinateType stdev = static_cast<LmaIsotopeFitter1D::Data*> (params) ->isotopes_stdev; //gsl_vector_get( x, 1 );
        CoordinateType mono_mz = gsl_vector_get( x, 1 );

        // iterate over all points of the signal
        for ( Size i = 0; i < n; ++i )
        {
            CoordinateType m = set[i].getPos();

            CoordinateType term1 = sqrt(2*Constants::PI)*stdev;
            CoordinateType termSum1 = 0.0;
            CoordinateType termSum2 = 0.0;
          //  CoordinateType termSum3 = 0.0;
            CoordinateType termExp = 0.0;

            for (Size j=0; j<isotopes_exact.size(); ++j)
            {
              termExp = exp(-pow(m-mono_mz-j*isotope_distance,2)/(2*stdev*stdev));
              termSum1 += isotopes_exact[j] * termExp;
              termSum2 += isotopes_exact[j] * termExp * ((m-mono_mz-j*isotope_distance)/(stdev*stdev));
           //   termSum3 += isotopes_exact[j] * termExp * (-1/stdev + (pow(m-mono_mz-j*isotope_distance,2)/(stdev*stdev*stdev)));
            }

           	CoordinateType f_a = (1/term1 * termSum1);
         //  	CoordinateType f_stdev = (A/term1 * termSum3);
           	CoordinateType f_mono_mz = (A/term1 * termSum2);

            // set the jacobian matrix
            gsl_matrix_set( J, i, 0, f_a/sigma );
       	 //   gsl_matrix_set( J, i, 1, f_stdev );
          	gsl_matrix_set( J, i, 1, f_mono_mz/sigma );
        }

        return GSL_SUCCESS;
    }

    Int LmaIsotopeFitter1D::evaluate_(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
    {
        LmaIsotopeFitter1D::residual_( x, params, f );
        LmaIsotopeFitter1D::jacobian_( x, params, J );

        return GSL_SUCCESS;
    }

    void LmaIsotopeFitter1D::printState_(Int iter, gsl_multifit_fdfsolver * s)
    {
        printf ( "iter: %4u x = % 15.8f % 15.8f |f(x)| = %g\n", iter,
                 gsl_vector_get( s->x, 0 ),
                 gsl_vector_get( s->x, 1 ),
                 //gsl_vector_get( s->x, 2 ),
                 gsl_blas_dnrm2( s->f ) );
    }

    void LmaIsotopeFitter1D::setInitialParameters_()
    {
        // compute the relative abundance of i-th isotopic peak
        isotopes_exact_.clear();
        typedef std::vector < DoubleReal > ContainerType;
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

        EmpiricalFormula formula(form);
        typedef IsotopeDistribution::iterator IsoIter;
        IsotopeDistribution isotope_distribution = formula.getIsotopeDistribution(max_isotope_);
        isotope_distribution.trimRight(trim_right_cutoff_);
        isotope_distribution.renormalize();

        // compute relative abundance of i-th isotopic peak of a peptide
        CoordinateType isotopes_mean = 0;
        Int i=0;
        for (	IsoIter iter = isotope_distribution.begin(); iter != isotope_distribution.end(); ++iter, ++i)
        {
          isotopes_exact_.push_back(iter->second);
          isotopes_mean += iter->second*i;
        }
        isotopes_mean *= isotope_distance_ / charge_;

        // compute monoisotopic mass
        if (monoisotopic_mz_ == 0.0)
        {
          monoisotopic_mz_ = mean_-isotopes_mean;
        }
    }

    LmaIsotopeFitter1D::QualityType LmaIsotopeFitter1D::fit1d(const RawDataArrayType& set, InterpolationModel*& model)
    {
        // Calculate bounding box
        min_ = max_ = set[0].getPos();
        for ( Size pos=1; pos < set.size(); ++pos)
        {
            CoordinateType tmp = set[pos].getPos();
            if ( min_ > tmp ) min_ = tmp;
            if ( max_ < tmp ) max_ = tmp;
        }

        if (monoisotopic_mz_ != 0.0) monoisotopic_mass_known_ = true;
        else monoisotopic_mass_known_ = false;

        // Enlarge the bounding box by a few multiples of the standard deviation
        {
          stdev1_ = sqrt ( statistics_.variance() ) * tolerance_stdev_box_;
          min_ -= stdev1_;
          max_ += stdev1_;
        }

        // Compute start parameters
        setInitialParameters_();

        // Set advanced parameters for residual_  und jacobian_ method
        LmaIsotopeFitter1D::Data d;
        d.n = set.size();
        d.set = set;
        d.isotopes_exact = isotopes_exact_;
        d.isotope_distance = isotope_distance_;
        // d.mono_known = monoisotopic_mass_known_;
        // d.monoisotopic_mz = monoisotopic_mz_;
        d.isotopes_stdev = isotope_stdev_;
        d.sigma = 0.1; // gaussian noise (standard deviation = 0.1)

        // The various _f, _df and _fdf of the LM loop return GSL_EDOM if any parameter is < 0

     	  // Optimize parameters with Levenberg-Marquardt algorithm (GLS)
        CoordinateType x_init[ 2 ] = { total_intensity_/100 , /*isotope_stdev_,*/ monoisotopic_mz_ };
        optimize_(set, 2, x_init, &(residual_), &(jacobian_), &(evaluate_), &d);

        // Set optimized parameters
        total_intensity_ = x_init[0];
        // isotope_stdev_ = x_init[1];
        monoisotopic_mz_ = x_init[1];

#ifdef DEBUG_FEATUREFINDER
        if ( getGslStatus_() != "success" )
        {
          std::cout << "status: " << getGslStatus_() << std::endl;
        }
#endif
        // build model
        if (charge_==0)
        {
          model = static_cast<InterpolationModel*> (Factory<BaseModel<1> >::create("GaussModel"));
          model->setInterpolationStep( interpolation_step_ );

          Param tmp;
          tmp.setValue( "bounding_box:min", min_ );
          tmp.setValue( "bounding_box:max", max_ );
          tmp.setValue( "statistics:variance", statistics_.variance() );
          tmp.setValue( "statistics:mean", statistics_.mean() );
          model->setParameters( tmp );
        }
        else
        {
          model = static_cast<InterpolationModel*> (Factory<BaseModel<1> >::create("LmaIsotopeModel"));

          Param iso_param = this->param_.copy( "isotope_model:", true );
          iso_param.removeAll("stdev");
          model->setParameters( iso_param );
          model->setInterpolationStep( interpolation_step_ );

          Param tmp;
          tmp.setValue( "total_intensity", total_intensity_ );
          tmp.setValue( "monoisotopic_mz", monoisotopic_mz_ );
          tmp.setValue( "bounding_box:min", min_ );
          tmp.setValue( "bounding_box:max", max_ );
          tmp.setValue( "statistics:mean", statistics_.mean() );
          tmp.setValue( "charge", static_cast<Int>( charge_ ) );
          tmp.setValue( "isotope:stdev", isotope_stdev_ );
		      tmp.setValue( "isotope:maximum", max_isotope_ );

          model->setParameters( tmp );
        }

  			// calculate pearson correlation
     		std::vector<Real> real_data;
        real_data.reserve(set.size());
        std::vector<Real> model_data;
        model_data.reserve(set.size());

        for (Size i=0; i < set.size(); ++i)
        {
           real_data.push_back(set[i].getIntensity());
           model_data.push_back( model->getIntensity( DPosition<1>(set[i].getPosition()) ) );
        }

        QualityType correlation = Math::pearsonCorrelationCoefficient(real_data.begin(), real_data.end(), model_data.begin(), model_data.end());
        if (boost::math::isnan(correlation)) correlation = -1.0;

        return correlation;
    }

    void LmaIsotopeFitter1D::updateMembers_()
    {
      LevMarqFitter1D::updateMembers_();

      total_intensity_ = param_.getValue("total_intensity");
      monoisotopic_mz_ = param_.getValue("monoisotopic_mass");

      statistics_.setVariance(param_.getValue("statistics:variance"));
      charge_ = param_.getValue("charge");
      isotope_stdev_ = param_.getValue("isotope:stdev");

      charge_ = param_.getValue("charge");
      mean_ = param_.getValue("statistics:mean");
      max_isotope_ = param_.getValue("isotope:maximum");

      trim_right_cutoff_ = param_.getValue("isotope:trim_right_cutoff");
      isotope_distance_ = param_.getValue("isotope:distance");

      averagine_[C] = param_.getValue("averagines:C");
      averagine_[H] = param_.getValue("averagines:H");
      averagine_[N] = param_.getValue("averagines:N");
      averagine_[O] = param_.getValue("averagines:O");
      averagine_[S] = param_.getValue("averagines:S");
    }

}
