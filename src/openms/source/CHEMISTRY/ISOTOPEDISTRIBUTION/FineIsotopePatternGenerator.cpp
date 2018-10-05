// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Hannes Rost $
// $Authors: Hannes Rost $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsoSpec.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/Element.h>

namespace OpenMS
{

  IsotopeDistribution FineIsotopePatternGenerator::run(const EmpiricalFormula& formula) const
  {
    IsoSpec algorithm(threshold_, absolute_);
#if 0
    // Use IsoSpec's isotopic tables
    algorithm.run(formula.toString());
#else
    // Use our own isotopic tables
    std::vector<int> isotopeNumbers, atomCounts;
    std::vector<std::vector<double> > isotopeMasses, isotopeProbabilities;

    // Iterate through all elements in the molecular formula
    for (auto elem : formula)
    {
      atomCounts.push_back(elem.second);

      std::vector<double> masses;
      std::vector<double> probs;
      for (auto iso : elem.first->getIsotopeDistribution())
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

    algorithm.run(isotopeNumbers, atomCounts, isotopeMasses, isotopeProbabilities);
#endif

    // Store the data in a IsotopeDistribution
    std::vector<Peak1D> c;
    c.reserve( algorithm.getMasses().size() );
    auto mit = algorithm.getMasses().cbegin();
    auto pit = algorithm.getProbabilities().cbegin();
    while (mit != algorithm.getMasses().cend())
    {
      c.emplace_back( Peak1D(*mit, *pit) );
      mit++; pit++;
    }

    IsotopeDistribution result;
    result.set(std::move(c));
    result.sortByMass();
    return result;
  }

}
