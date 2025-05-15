// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/ExtendedIsotopeModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

namespace OpenMS
{
  ExtendedIsotopeModel::ExtendedIsotopeModel() :
    InterpolationModel(),
    charge_(0),
    monoisotopic_mz_(0.0)
  {
    setName("ExtendedIsotopeModel");

    defaults_.setValue("averagines:C", 0.04443989f, "Number of C atoms per Dalton of mass.", {"advanced"});
    defaults_.setValue("averagines:H", 0.06981572f, "Number of H atoms per Dalton of mass.", {"advanced"});
    defaults_.setValue("averagines:N", 0.01221773f, "Number of N atoms per Dalton of mass.", {"advanced"});
    defaults_.setValue("averagines:O", 0.01329399f, "Number of O atoms per Dalton of mass.", {"advanced"});
    defaults_.setValue("averagines:S", 0.00037525f, "Number of S atoms per Dalton of mass.", {"advanced"});

    defaults_.setValue("isotope:trim_right_cutoff", 0.001, "Cutoff in averagine distribution, trailing isotopes below this relative intensity are not considered.", {"advanced"});
    defaults_.setValue("isotope:maximum", 100, "Maximum isotopic rank to be considered.", {"advanced"});
    defaults_.setValue("isotope:distance", 1.000495, "Distance between consecutive isotopic peaks.", {"advanced"});
    defaults_.setValue("isotope:stdev", 0.1, "Standard deviation of gaussian applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.", {"advanced"});
    defaults_.setValue("charge", 1, "Charge state of the model.", {"advanced"});
    defaults_.setValue("isotope:monoisotopic_mz", 1.0, "Monoisotopic m/z of the model.", {"advanced"});

    defaultsToParam_();
  }

  ExtendedIsotopeModel::ExtendedIsotopeModel(const ExtendedIsotopeModel & source) :
    InterpolationModel(source)
  {
    setParameters(source.getParameters());
    updateMembers_();
  }

  ExtendedIsotopeModel::~ExtendedIsotopeModel() = default;

  ExtendedIsotopeModel & ExtendedIsotopeModel::operator=(const ExtendedIsotopeModel & source)
  {
    if (&source == this)
    {
      return *this;
    }
    InterpolationModel::operator=(source);
    setParameters(source.getParameters());
    updateMembers_();

    return *this;
  }

  void ExtendedIsotopeModel::setSamples()
  {
    // MAGIC alert, num stdev for smooth table for normal distribution
    CoordinateType normal_widening_num_stdev = 4.;
    // Actual width for values in the smooth table for normal distribution
    CoordinateType normal_widening_width = isotope_stdev_ * normal_widening_num_stdev;

    typedef std::vector<double> ContainerType;
    ContainerType isotopes_exact;
    CoordinateType mass = monoisotopic_mz_ * charge_;

    Int C_num = Int(0.5 + mass * averagine_[C]);
    Int N_num = Int(0.5 + mass * averagine_[N]);
    Int O_num = Int(0.5 + mass * averagine_[O]);
    Int H_num = Int(0.5 + mass * averagine_[H]);
    Int S_num = Int(0.5 + mass * averagine_[S]);

    String form("");
    if (C_num)
    {
      form.append("C").append(String(C_num));
    }
    if (H_num)
    {
      form.append("H").append(String(H_num));
    }
    if (N_num)
    {
      form.append("N").append(String(N_num));
    }
    if (O_num)
    {
      form.append("O").append(String(O_num));
    }
    if (S_num)
    {
      form.append("S").append(String(S_num));
    }
    EmpiricalFormula formula(form);
    IsotopeDistribution isotope_distribution = formula.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
    isotope_distribution.trimRight(trim_right_cutoff_);
    isotope_distribution.renormalize();

    // compute the average mass (-offset)
    for (const Peak1D& peak : isotope_distribution)
    {
      isotopes_exact.push_back(peak.getIntensity());
    }

    // "stretch" the averagine isotope distribution
    Size isotopes_exact_size = isotopes_exact.size();
    isotopes_exact.resize(Size((isotopes_exact_size - 1)
                               * isotope_distance_ / interpolation_step_ + 1.6));                             // round up a bit more

    for (Size i = isotopes_exact_size - 1; i; --i)
    {
      // we don't need to move the 0-th entry
      isotopes_exact[Size(CoordinateType(i) *
                          isotope_distance_ / interpolation_step_ / charge_ + 0.5)]
        =   isotopes_exact[i];
      isotopes_exact[i] = 0;
    }

    // compute the normal distribution (to be added for widening the averagine isotope distribution)
    Math::BasicStatistics<> normal_widening_model;
    normal_widening_model.setSum(1);
    normal_widening_model.setMean(0);
    normal_widening_model.setVariance(isotope_stdev_ * isotope_stdev_);
    // fill a container with CoordinateType points
    ContainerType normal_widening_coordinate;
    for (double coord = -normal_widening_width;
         coord <= normal_widening_width;
         coord += interpolation_step_
         )
    {
      normal_widening_coordinate.push_back(coord);
    }
    // compute normal approximation at these CoordinateType points
    ContainerType normal_widening;
    normal_widening_model.normalApproximation(normal_widening, normal_widening_coordinate);

    // fill linear interpolation
    const ContainerType & left = isotopes_exact;
    const ContainerType & right = normal_widening;
    ContainerType & result = interpolation_.getData();
    result.clear();

    Int rMax = std::min(Int(left.size() + right.size() - 1), Int(2 * normal_widening_width / interpolation_step_ * max_isotope_ + 1));
    result.resize(rMax, 0);

    // we loop backwards because then the small products tend to come first
    // (for better numerics)
    for (SignedSize i = left.size() - 1; i >= 0; --i)
    {
      if (left[i] == 0)
      {
        continue;
      }
      for (SignedSize j = std::min<SignedSize>(rMax - i, right.size()) - 1; j >= 0; --j)
      {
        result[i + j] += left[i] * right[j];
      }
    }

    // set interpolation
    interpolation_.setMapping(interpolation_step_, normal_widening_width / interpolation_step_, monoisotopic_mz_);

    // scale data so that integral over distribution equals one
    // multiply sum by interpolation_step_ -> rectangular approximation of integral
    IntensityType factor = scaling_ / interpolation_step_ /
                           std::accumulate(result.begin(), result.end(), IntensityType(0));

    for (ContainerType::iterator iter = result.begin(); iter != result.end(); ++iter)
    {
      *iter *= factor;
    }

  }

  void ExtendedIsotopeModel::setOffset(CoordinateType offset)
  {
    double diff = offset - getInterpolation().getOffset();
    monoisotopic_mz_ += diff;

    InterpolationModel::setOffset(offset);

    param_.setValue("isotope:monoisotopic_mz", monoisotopic_mz_);
  }

  ExtendedIsotopeModel::CoordinateType ExtendedIsotopeModel::getOffset()
  {
    return getInterpolation().getOffset();
  }

  UInt ExtendedIsotopeModel::getCharge() const
  {
    return charge_;
  }

  ExtendedIsotopeModel::CoordinateType ExtendedIsotopeModel::getCenter() const
  {
    return monoisotopic_mz_;
  }

  void ExtendedIsotopeModel::updateMembers_()
  {
    InterpolationModel::updateMembers_();

    charge_ = param_.getValue("charge");
    isotope_stdev_ = param_.getValue("isotope:stdev");
    monoisotopic_mz_ = param_.getValue("isotope:monoisotopic_mz");
    max_isotope_ = param_.getValue("isotope:maximum");
    trim_right_cutoff_ = param_.getValue("isotope:trim_right_cutoff");
    isotope_distance_ = param_.getValue("isotope:distance");

    averagine_[C] = param_.getValue("averagines:C");
    averagine_[H] = param_.getValue("averagines:H");
    averagine_[N] = param_.getValue("averagines:N");
    averagine_[O] = param_.getValue("averagines:O");
    averagine_[S] = param_.getValue("averagines:S");

    setSamples();
  }

}
