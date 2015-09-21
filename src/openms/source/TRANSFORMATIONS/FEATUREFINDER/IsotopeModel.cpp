// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#include <boost/math/distributions/cauchy.hpp>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <numeric>

namespace OpenMS
{
  IsotopeModel::IsotopeModel() :
    InterpolationModel(),
    charge_(0),
    monoisotopic_mz_(0.0)
  {
    setName(getProductName());

    defaults_.setValue("averagines:C", 0.04443989f, "Number of C atoms per Dalton of mass.", ListUtils::create<String>("advanced"));
    defaults_.setValue("averagines:H", 0.06981572f, "Number of H atoms per Dalton of mass.", ListUtils::create<String>("advanced"));
    defaults_.setValue("averagines:N", 0.01221773f, "Number of N atoms per Dalton of mass.", ListUtils::create<String>("advanced"));
    defaults_.setValue("averagines:O", 0.01329399f, "Number of O atoms per Dalton of mass.", ListUtils::create<String>("advanced"));
    defaults_.setValue("averagines:S", 0.00037525f, "Number of S atoms per Dalton of mass.", ListUtils::create<String>("advanced"));
    defaults_.setValue("isotope:trim_right_cutoff", 0.001, "Cutoff in averagine distribution, trailing isotopes below this relative intensity are not considered.", ListUtils::create<String>("advanced"));
    defaults_.setValue("isotope:maximum", 100, "Maximum isotopic rank to be considered.", ListUtils::create<String>("advanced"));
    defaults_.setValue("isotope:distance", 1.000495, "Distance between consecutive isotopic peaks.", ListUtils::create<String>("advanced"));


    defaults_.setValue("isotope:mode:mode", "Gaussian", "Peak Shape used around each isotope peak.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("isotope:mode:mode", ListUtils::create<String>("Gaussian,Lorentzian"));
    defaults_.setValue("isotope:mode:LorentzFWHM", 0.3, "Full width of the Lorentzian (Cauchy) function applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.", ListUtils::create<String>("advanced"));
    defaults_.setValue("isotope:mode:GaussianSD", 0.1, "Standard deviation of Gaussian applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.", ListUtils::create<String>("advanced"));


    defaults_.setValue("charge", 1, "Charge state of the model.", ListUtils::create<String>("advanced"));
    defaults_.setValue("statistics:mean", 0.0, "Centroid m/z (as opposed to monoisotopic m/z).", ListUtils::create<String>("advanced"));

    defaultsToParam_();
  }

  IsotopeModel::IsotopeModel(const IsotopeModel & source) :
    InterpolationModel(source)
  {
    setParameters(source.getParameters());
    updateMembers_();
  }

  IsotopeModel::~IsotopeModel()
  {
  }

  IsotopeModel & IsotopeModel::operator=(const IsotopeModel & source)
  {
    if (&source == this)
      return *this;

    InterpolationModel::operator=(source);
    setParameters(source.getParameters());
    updateMembers_();

    return *this;
  }

  EmpiricalFormula IsotopeModel::getFormula()
  {
    CoordinateType mass = mean_ * charge_;

    Int C_num = Int(0.5 + mass * averagine_[C]);
    Int N_num = Int(0.5 + mass * averagine_[N]);
    Int O_num = Int(0.5 + mass * averagine_[O]);
    Int H_num = Int(0.5 + mass * averagine_[H]);
    Int S_num = Int(0.5 + mass * averagine_[S]);

    String form;
    if (C_num)
      form.append("C").append(String(C_num));
    if (H_num)
      form.append("H").append(String(H_num));
    if (N_num)
      form.append("N").append(String(N_num));
    if (O_num)
      form.append("O").append(String(O_num));
    if (S_num)
      form.append("S").append(String(S_num));

    return EmpiricalFormula(form);
  }

  const IsotopeDistribution & IsotopeModel::getIsotopeDistribution() const
  {
    return isotope_distribution_;
  }

  void IsotopeModel::setSamples(const EmpiricalFormula & formula)
  {
    typedef std::vector<double> ContainerType;
    ContainerType isotopes_exact;

    isotope_distribution_ = formula.getIsotopeDistribution(max_isotope_);

    isotope_distribution_.trimRight(trim_right_cutoff_);
    isotope_distribution_.renormalize();

    // compute the average mass (-offset)
    CoordinateType isotopes_mean = 0;
    {
      Int cnt = 0;
      for (IsotopeDistribution::iterator iter = isotope_distribution_.begin();
           iter != isotope_distribution_.end(); ++iter, ++cnt)
      {
        isotopes_exact.push_back(iter->second);
        isotopes_mean += iter->second * cnt;
      }
      isotopes_mean *= isotope_distance_ / charge_;
    }
    // (Need not divide by sum of probabilities, which is 1.)

    ///
    // "stretch" the averagine isotope distribution (so we can add datapoints between isotope peaks)
    ///
    size_t isotopes_exact_size = isotopes_exact.size();
    isotopes_exact.resize(size_t((isotopes_exact_size - 1) * isotope_distance_ / interpolation_step_ + 1.6)); // round up a bit more

    for (Size i = isotopes_exact_size - 1; i; --i)
    {
      // we don't need to move the 0-th entry
      isotopes_exact[size_t(CoordinateType(i) * isotope_distance_ / interpolation_step_ / charge_ + 0.5)]
        =   isotopes_exact[i];
      isotopes_exact[i] = 0;
    }

    ////
    // compute the Gaussian/Cauchy distribution (to be added for widening the averagine isotope distribution)
    ////
    ContainerType peak_shape_values_y;
    // fill a container with CoordinateType points (x values)
    CoordinateType peak_width = 0.0;
    if (param_.getValue("isotope:mode:mode") == "Gaussian")
    {
      // Actual width for values in the smooth table for normal distribution
      peak_width = isotope_stdev_ * 4.0;  // MAGIC alert, num stdev for smooth table for normal distribution
      ContainerType peak_shape_values_x;
      for (double coord = -peak_width; coord <= peak_width;
           coord += interpolation_step_)
      {
        peak_shape_values_x.push_back(coord);
      }
      // compute normal approximation at these CoordinateType points (y values)
      Math::BasicStatistics<> normal_widening_model;
      normal_widening_model.setSum(1);
      normal_widening_model.setMean(0);
      normal_widening_model.setVariance(isotope_stdev_ * isotope_stdev_);
      normal_widening_model.normalApproximation(peak_shape_values_y, peak_shape_values_x);
    }
    else if (param_.getValue("isotope:mode:mode") == "Lorentzian")
    {
      peak_width = isotope_lorentz_fwhm_ * 8.0; // MAGIC alert: Lorentzian has infinite support, but we need to stop sampling at some point: 8*FWHM
      for (double coord = -peak_width; coord <= peak_width;
           coord += interpolation_step_)
      {
        boost::math::cauchy_distribution<double> cauchy(0., isotope_lorentz_fwhm_ / 2.0);
        double x = boost::math::pdf(cauchy, coord);
        peak_shape_values_y.push_back(x); //cauchy is using HWHM not FWHM
      }
    }

    ///
    // fold the Gaussian/Lorentzian at each averagine peak, i.e. fill linear interpolation
    ///
    const ContainerType & left = isotopes_exact;
    const ContainerType & right = peak_shape_values_y;
    ContainerType & result = interpolation_.getData();
    result.clear();

    SignedSize r_max = std::min(SignedSize(left.size() + right.size() - 1),
                                SignedSize(2 * peak_width / interpolation_step_ * max_isotope_ + 1));
    result.resize(r_max, 0);

    // we loop backwards because then the small products tend to come first
    // (for better numerics)
    for (SignedSize i = left.size() - 1; i >= 0; --i)
    {
      if (left[i] == 0)
        continue;
      for (SignedSize j = std::min(r_max - i, SignedSize(right.size())) - 1; j >= 0; --j)
      {
        result[i + j] += left[i] * right[j];
      }
    }

    monoisotopic_mz_ = mean_ - isotopes_mean;
    interpolation_.setMapping(interpolation_step_, peak_width / interpolation_step_, monoisotopic_mz_);

    //std::cerr << "mono now: " << monoisotopic_mz_ << " mono easy: " << formula.getMonoWeight()/formula.getCharge() << "\n";

    // scale data so that integral over distribution equals one
    // multiply sum by interpolation_step_ -> rectangular approximation of integral
    IntensityType factor = scaling_ / (interpolation_step_ * std::accumulate(result.begin(), result.end(), IntensityType(0)));
    for (ContainerType::iterator iter = result.begin(); iter != result.end(); ++iter)
    {
      *iter *= factor;
    }
  }

  void IsotopeModel::setOffset(CoordinateType offset)
  {
    double diff = offset - getInterpolation().getOffset();
    mean_ += diff;
    monoisotopic_mz_ += diff;

    InterpolationModel::setOffset(offset);

    param_.setValue("statistics:mean", mean_);
  }

  IsotopeModel::CoordinateType IsotopeModel::getOffset()
  {
    return getInterpolation().getOffset();
  }

  UInt IsotopeModel::getCharge()
  {
    return charge_;
  }

  IsotopeModel::CoordinateType IsotopeModel::getCenter() const
  {
    return monoisotopic_mz_;
  }

  void IsotopeModel::updateMembers_()
  {
    InterpolationModel::updateMembers_();

    charge_ = param_.getValue("charge");
    isotope_stdev_ = param_.getValue("isotope:mode:GaussianSD");
    isotope_lorentz_fwhm_ = param_.getValue("isotope:mode:LorentzFWHM");
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
