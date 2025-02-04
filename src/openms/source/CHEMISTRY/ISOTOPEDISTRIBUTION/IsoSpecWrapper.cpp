// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Hannes Rost, Michał Startek, Mateusz Łącki $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsoSpecWrapper.h>

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <iterator>
#include <memory>
#include <string>
#include <utility>

// Override IsoSpec's use of mmap whenever it is available
#define ISOSPEC_GOT_SYSTEM_MMAN false
#define ISOSPEC_GOT_MMAN false
#define ISOSPEC_BUILDING_OPENMS true

// TODO Fix this weird way of including. Just make a library out of it,
//  link to it and declare it PUBLIC such that it gets linked to dependents of OpenMS lib
// But since it is PUBLIC, you should export the library also for installation (see evergreen thirdparty)
#include "IsoSpec/allocator.cpp"
#include "IsoSpec/dirtyAllocator.cpp"
#include "IsoSpec/isoSpec++.cpp"
#include "IsoSpec/isoMath.cpp"
#include "IsoSpec/marginalTrek++.cpp"
#include "IsoSpec/operators.cpp"
#include "IsoSpec/element_tables.cpp"
#include "IsoSpec/misc.cpp"
#include "IsoSpec/fasta.cpp"

using namespace std;
using namespace IsoSpec;

namespace OpenMS
{
  /// Convert an set of isotope probabiities to IsoSpec input
  Iso _OMS_IsoFromParameters(const std::vector<int>& isotopeNr,
                             const std::vector<int>& atomCounts,
                             const std::vector<std::vector<double> >& isotopeMasses,
                             const std::vector<std::vector<double> >& isotopeProbabilities)
  {
    OPENMS_PRECONDITION(isotopeNr.size() == atomCounts.size(), "Vectors need to be of the same size")
    OPENMS_PRECONDITION(isotopeNr.size() == isotopeMasses.size(), "Vectors need to be of the same size")
    OPENMS_PRECONDITION(isotopeNr.size() == isotopeProbabilities.size(), "Vectors need to be of the same size")

    // Check that all probabilities are non-zero
    if (!std::all_of(std::begin(isotopeProbabilities), std::end(isotopeProbabilities), [](std::vector<double> prob){ 
            return std::all_of(std::begin(prob), std::end(prob), [](double p){return p > 0.0;});
          })) 
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       std::string("All probabilities need to be larger than zero").c_str());
    }

    int dimNumber = isotopeNr.size();

    // Convert vector of vector to double**
    std::unique_ptr<const double*[]> IM(new const double*[dimNumber]);
    std::unique_ptr<const double*[]> IP(new const double*[dimNumber]);
    for (int i = 0; i < dimNumber; i++)
    {
      IM[i] = isotopeMasses[i].data();
      IP[i] = isotopeProbabilities[i].data();
    }

    // IsoSpec will *copy* these values, so once provided we are safe to
    // destroy them here on the OpenMS side
    Iso ret(dimNumber, isotopeNr.data(), atomCounts.data(), IM.get(), IP.get());
    
    return ret;
  }

  /// Convert an OpenMS EmpiricalFormula to the input format for IsoSpec
  Iso _OMS_IsoFromEmpiricalFormula(const EmpiricalFormula& formula) 
  {
    // Use our own isotopic tables
    std::vector<int> isotopeNumbers, atomCounts;
    std::vector<std::vector<double> > isotopeMasses, isotopeProbabilities;

    // Iterate through all elements in the molecular formula
    for (const auto& elem : formula)
    {
      atomCounts.push_back(elem.second);

      std::vector<double> masses;
      std::vector<double> probs;
      for (const auto& iso : elem.first->getIsotopeDistribution())
      {
        if (iso.getIntensity() <= 0.0) continue; // Note: there will be a segfault if one of the intensities is zero!
        masses.push_back(iso.getMZ());
        probs.push_back(iso.getIntensity());
      }

      // For each element store how many isotopes it has and their masses/probabilities
      isotopeNumbers.push_back( masses.size() );
      isotopeMasses.push_back(masses);
      isotopeProbabilities.push_back(probs);
    }

    // Store the data in a IsotopeDistribution
    return _OMS_IsoFromParameters(isotopeNumbers, atomCounts, isotopeMasses, isotopeProbabilities);
  }

  IsoSpecThresholdGeneratorWrapper::IsoSpecThresholdGeneratorWrapper(const std::vector<int>& isotopeNr,
                    const std::vector<int>& atomCounts,
                    const std::vector<std::vector<double> >& isotopeMasses,
                    const std::vector<std::vector<double> >& isotopeProbabilities,
                    double threshold,
                    bool absolute) :
  ITG(std::make_unique<IsoSpec::IsoThresholdGenerator>(
    _OMS_IsoFromParameters(isotopeNr, atomCounts, isotopeMasses, isotopeProbabilities), 
    threshold, 
    absolute))
  {};

  IsoSpecThresholdGeneratorWrapper::IsoSpecThresholdGeneratorWrapper(const EmpiricalFormula& formula,
                    double threshold,
                    bool absolute) :
  ITG(std::make_unique<IsoSpec::IsoThresholdGenerator>(_OMS_IsoFromEmpiricalFormula(formula), threshold, absolute))
  {};

  bool IsoSpecThresholdGeneratorWrapper::nextConf() { return ITG->advanceToNextConfiguration(); };
  Peak1D IsoSpecThresholdGeneratorWrapper::getConf() { return Peak1D(ITG->mass(), ITG->prob()); };
  double IsoSpecThresholdGeneratorWrapper::getMass() { return ITG->mass(); };
  double IsoSpecThresholdGeneratorWrapper::getIntensity() { return ITG->prob(); };
  double IsoSpecThresholdGeneratorWrapper::getLogIntensity() { return ITG->lprob(); };

  // in this special case it needs to go in cpp file (see e.g., https://stackoverflow.com/questions/38242200/where-should-a-default-destructor-c11-style-go-header-or-cpp)
  IsoSpecThresholdGeneratorWrapper::~IsoSpecThresholdGeneratorWrapper() = default; 

//  --------------------------------------------------------------------------------

  IsoSpecTotalProbGeneratorWrapper::IsoSpecTotalProbGeneratorWrapper(const std::vector<int>& isotopeNr,
                    const std::vector<int>& atomCounts,
                    const std::vector<std::vector<double> >& isotopeMasses,
                    const std::vector<std::vector<double> >& isotopeProbabilities,
                    double total_prob_hint) :
  ILG(std::make_unique<IsoSpec::IsoLayeredGenerator>(_OMS_IsoFromParameters(isotopeNr, atomCounts, isotopeMasses, isotopeProbabilities), 1024, 1024, true, total_prob_hint))
  {};

  IsoSpecTotalProbGeneratorWrapper::IsoSpecTotalProbGeneratorWrapper(const EmpiricalFormula& formula,
                    double total_prob_hint) :
  ILG(std::make_unique<IsoSpec::IsoLayeredGenerator>(_OMS_IsoFromEmpiricalFormula(formula), 1024, 1024, true, total_prob_hint))
  {};

  IsoSpecTotalProbGeneratorWrapper::~IsoSpecTotalProbGeneratorWrapper() = default;

  bool IsoSpecTotalProbGeneratorWrapper::nextConf() { return ILG->advanceToNextConfiguration(); };
  Peak1D IsoSpecTotalProbGeneratorWrapper::getConf() { return Peak1D(ILG->mass(), ILG->prob()); };
  double IsoSpecTotalProbGeneratorWrapper::getMass() { return ILG->mass(); };
  double IsoSpecTotalProbGeneratorWrapper::getIntensity() { return ILG->prob(); };
  double IsoSpecTotalProbGeneratorWrapper::getLogIntensity() { return ILG->lprob(); };

//  --------------------------------------------------------------------------------

  IsoSpecOrderedGeneratorWrapper::IsoSpecOrderedGeneratorWrapper(const std::vector<int>& isotopeNr,
                    const std::vector<int>& atomCounts,
                    const std::vector<std::vector<double> >& isotopeMasses,
                    const std::vector<std::vector<double> >& isotopeProbabilities) :
   IOG(std::make_unique<IsoSpec::IsoOrderedGenerator>(_OMS_IsoFromParameters(isotopeNr, atomCounts, isotopeMasses, isotopeProbabilities)))
  {};

  IsoSpecOrderedGeneratorWrapper::IsoSpecOrderedGeneratorWrapper(const EmpiricalFormula& formula) :
    IOG(std::make_unique<IsoSpec::IsoOrderedGenerator>(_OMS_IsoFromEmpiricalFormula(formula)))
  {};

  IsoSpecOrderedGeneratorWrapper::~IsoSpecOrderedGeneratorWrapper() = default; // needs to be in cpp file because of incomplete types!

  bool IsoSpecOrderedGeneratorWrapper::nextConf() { return IOG->advanceToNextConfiguration(); };
  Peak1D IsoSpecOrderedGeneratorWrapper::getConf() { return Peak1D(IOG->mass(), IOG->prob()); };
  double IsoSpecOrderedGeneratorWrapper::getMass() { return IOG->mass(); };
  double IsoSpecOrderedGeneratorWrapper::getIntensity() { return IOG->prob(); };
  double IsoSpecOrderedGeneratorWrapper::getLogIntensity() { return IOG->lprob(); };

//  --------------------------------------------------------------------------------

  IsoSpecThresholdWrapper::IsoSpecThresholdWrapper(const std::vector<int>& isotopeNr,
                    const std::vector<int>& atomCounts,
                    const std::vector<std::vector<double> >& isotopeMasses,
                    const std::vector<std::vector<double> >& isotopeProbabilities,
                    double threshold,
                    bool absolute) :
  ITG(std::make_unique<IsoSpec::IsoThresholdGenerator>(
    _OMS_IsoFromParameters(isotopeNr, atomCounts, isotopeMasses, isotopeProbabilities), 
    threshold, 
    absolute))
  {}

  IsoSpecThresholdWrapper::IsoSpecThresholdWrapper(const EmpiricalFormula& formula,
                    double threshold,
                    bool absolute) :                    
  ITG(std::make_unique<IsoSpec::IsoThresholdGenerator>(
    _OMS_IsoFromEmpiricalFormula(formula), 
    threshold, 
    absolute))
  {};

  IsotopeDistribution IsoSpecThresholdWrapper::run()
  {
    std::vector<Peak1D> distribution;
    distribution.reserve(ITG->count_confs());

    ITG->reset();

    while (ITG->advanceToNextConfiguration())
        distribution.emplace_back(Peak1D(ITG->mass(), ITG->prob()));

    IsotopeDistribution ID;

    ID.set(std::move(distribution));

    return ID;
  }

  IsoSpecThresholdWrapper::~IsoSpecThresholdWrapper() = default;
  
//  --------------------------------------------------------------------------------


  IsoSpecTotalProbWrapper::IsoSpecTotalProbWrapper(const std::vector<int>& isotopeNr,
                    const std::vector<int>& atomCounts,
                    const std::vector<std::vector<double> >& isotopeMasses,
                    const std::vector<std::vector<double> >& isotopeProbabilities,
                    double _total_prob,
                    bool _do_p_trim) :
  ILG(std::make_unique<IsoSpec::IsoLayeredGenerator>(
    _OMS_IsoFromParameters(isotopeNr, atomCounts, isotopeMasses, isotopeProbabilities), 
    1024, 
    1024, 
    true, 
    _total_prob)),
  target_prob(_total_prob),
  do_p_trim(_do_p_trim)
  {};

  IsoSpecTotalProbWrapper::IsoSpecTotalProbWrapper(const EmpiricalFormula& formula,
                    double _total_prob,
                    bool _do_p_trim) :
  ILG(std::make_unique<IsoSpec::IsoLayeredGenerator>(_OMS_IsoFromEmpiricalFormula(formula), 1024, 1024, true, _total_prob)),
  target_prob(_total_prob),
  do_p_trim(_do_p_trim)
  {};

  IsoSpecTotalProbWrapper::~IsoSpecTotalProbWrapper() = default;

  IsotopeDistribution IsoSpecTotalProbWrapper::run()
  {
    std::vector<Peak1D> distribution;
    // There is no sensible way to precalculate the number of configurations 
    // in IsoLayeredGenerator

    double acc_prob = 0.0;

    while (acc_prob < target_prob && ILG->advanceToNextConfiguration())
    {
        double p = ILG->prob();
        acc_prob += p;
        distribution.emplace_back(Peak1D(ILG->mass(), p));
    }

    if (do_p_trim)
    {
        // the p_trim: extract the rest of the last layer, and perform quickselect

        while (ILG->advanceToNextConfigurationWithinLayer())
            distribution.emplace_back(Peak1D(ILG->mass(), ILG->prob()));

        size_t start = 0;
        size_t end = distribution.size();
        double sum_to_start = 0.0;

        while (start < end)
        {
            // Partition part
            size_t pivot = start + (end-start)/2; // middle
            double pprob = distribution[pivot].getIntensity();
            std::swap(distribution[pivot], distribution[end-1]);

            double new_csum = sum_to_start;

            size_t loweridx = start;
            for (size_t ii = start; ii < end-1; ii++)
                if (distribution[ii].getIntensity() > pprob)
                {
                    std::swap(distribution[ii], distribution[loweridx]);
                    new_csum += distribution[loweridx].getIntensity();
                    loweridx++;
                }

            std::swap(distribution[end-1], distribution[loweridx]);

            // Selection part
            if (new_csum < target_prob)
            {
                start = loweridx + 1;
                sum_to_start = new_csum + distribution[loweridx].getIntensity();
            }
            else
                end = loweridx;
        }
        distribution.resize(end);
    }

    IsotopeDistribution ID;
    ID.set(std::move(distribution));
    return ID;
  }
}  // namespace OpenMS
