// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#define ISOSPEC_GOT_SYSTEM_MMAN false
#define ISOSPEC_GOT_MMAN false
#define ISOSPEC_BUILDING_OPENMS true

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
    for (int i=0; i<dimNumber; i++)
    {
      IM[i] = isotopeMasses[i].data();
      IP[i] = isotopeProbabilities[i].data();
    }

    Iso ret(dimNumber, isotopeNr.data(), atomCounts.data(), IM.get(), IP.get());
    
    return ret;
  }

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
  ITG(_OMS_IsoFromParameters(isotopeNr, atomCounts, isotopeMasses, isotopeProbabilities), threshold, absolute)
  {};

  IsoSpecThresholdGeneratorWrapper::IsoSpecThresholdGeneratorWrapper(const EmpiricalFormula& formula,
                    double threshold,
                    bool absolute) :
  ITG(_OMS_IsoFromEmpiricalFormula(formula), threshold, absolute)
  {};


//  --------------------------------------------------------------------------------



  IsoSpecTotalProbGeneratorWrapper::IsoSpecTotalProbGeneratorWrapper(const std::vector<int>& isotopeNr,
                    const std::vector<int>& atomCounts,
                    const std::vector<std::vector<double> >& isotopeMasses,
                    const std::vector<std::vector<double> >& isotopeProbabilities,
                    double total_prob_hint) :
  ILG(_OMS_IsoFromParameters(isotopeNr, atomCounts, isotopeMasses, isotopeProbabilities), 1024, 1024, true, total_prob_hint)
  {};

  IsoSpecTotalProbGeneratorWrapper::IsoSpecTotalProbGeneratorWrapper(const EmpiricalFormula& formula,
                    double total_prob_hint) :
  ILG(_OMS_IsoFromEmpiricalFormula(formula), 1024, 1024, true, total_prob_hint)
  {};


//  --------------------------------------------------------------------------------


  IsoSpecOrderedGeneratorWrapper::IsoSpecOrderedGeneratorWrapper(const std::vector<int>& isotopeNr,
                    const std::vector<int>& atomCounts,
                    const std::vector<std::vector<double> >& isotopeMasses,
                    const std::vector<std::vector<double> >& isotopeProbabilities) :
  IOG(_OMS_IsoFromParameters(isotopeNr, atomCounts, isotopeMasses, isotopeProbabilities))
  {};

  IsoSpecOrderedGeneratorWrapper::IsoSpecOrderedGeneratorWrapper(const EmpiricalFormula& formula) :
  IOG(_OMS_IsoFromEmpiricalFormula(formula))
  {};

//  --------------------------------------------------------------------------------

  IsoSpecThresholdWrapper::IsoSpecThresholdWrapper(const std::vector<int>& isotopeNr,
                    const std::vector<int>& atomCounts,
                    const std::vector<std::vector<double> >& isotopeMasses,
                    const std::vector<std::vector<double> >& isotopeProbabilities,
                    double threshold,
                    bool absolute) :
  ITG(_OMS_IsoFromParameters(isotopeNr, atomCounts, isotopeMasses, isotopeProbabilities), threshold, absolute)
  {};

  IsoSpecThresholdWrapper::IsoSpecThresholdWrapper(const EmpiricalFormula& formula,
                    double threshold,
                    bool absolute) :
  ITG(_OMS_IsoFromEmpiricalFormula(formula), threshold, absolute)
  {};


  IsotopeDistribution IsoSpecThresholdWrapper::run()
  {
    std::vector<Peak1D> distribution;
    distribution.reserve(ITG.count_confs());

    ITG.reset();

    while (ITG.advanceToNextConfiguration())
        distribution.emplace_back(Peak1D(ITG.mass(), ITG.prob()));

    IsotopeDistribution ID;

    ID.set(std::move(distribution));

    return ID;
  }

//  --------------------------------------------------------------------------------


  IsoSpecTotalProbWrapper::IsoSpecTotalProbWrapper(const std::vector<int>& isotopeNr,
                    const std::vector<int>& atomCounts,
                    const std::vector<std::vector<double> >& isotopeMasses,
                    const std::vector<std::vector<double> >& isotopeProbabilities,
                    double _total_prob,
                    bool _do_p_trim) :
  ILG(_OMS_IsoFromParameters(isotopeNr, atomCounts, isotopeMasses, isotopeProbabilities), 1024, 1024, true, _total_prob),
  target_prob(_total_prob),
  do_p_trim(_do_p_trim)
  {};

  IsoSpecTotalProbWrapper::IsoSpecTotalProbWrapper(const EmpiricalFormula& formula,
                    double _total_prob,
                    bool _do_p_trim) :
  ILG(_OMS_IsoFromEmpiricalFormula(formula), 1024, 1024, true, _total_prob),
  target_prob(_total_prob),
  do_p_trim(_do_p_trim)
  {};


  IsotopeDistribution IsoSpecTotalProbWrapper::run()
  {
    std::vector<Peak1D> distribution;
    // There is no sensible way to precalculate the number of configurations 
    // in IsoLayeredGenerator

    double acc_prob = 0.0;

    while (acc_prob < target_prob && ILG.advanceToNextConfiguration())
    {
        double p = ILG.prob();
        acc_prob += p;
        distribution.emplace_back(Peak1D(ILG.mass(), p));
    }

    if (do_p_trim)
    {
        // the p_trim: extract the rest of the last layer, and perform quickselect

        while (ILG.advanceToNextConfigurationWithinLayer())
            distribution.emplace_back(Peak1D(ILG.mass(), ILG.prob()));

        size_t start = 0;
        size_t end = distribution.size();
        double sum_to_start = 0.0;

        while (start < end)
        {
            // Partition part
            size_t len = end - start;
            size_t pivot = len/2 + start;
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
