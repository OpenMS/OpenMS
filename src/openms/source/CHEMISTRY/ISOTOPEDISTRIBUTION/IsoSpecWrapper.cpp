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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsoSpecWrapper.h>

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <iterator>
#include <string>

// Override IsoSpec's use of mmap whenever it is available
#define ISOSPEC_TEST_GOT_SYSTEM_MMAN false
#define ISOSPEC_TEST_GOT_MMAN false

#include "IsoSpec/allocator.cpp"
#include "IsoSpec/dirtyAllocator.cpp"
#include "IsoSpec/isoSpec++.cpp"
#include "IsoSpec/isoMath.cpp"
#include "IsoSpec/marginalTrek++.cpp"
#include "IsoSpec/operators.cpp"
#include "IsoSpec/element_tables.cpp"
#include "IsoSpec/misc.cpp"

using namespace std;
using namespace IsoSpec;

namespace OpenMS
{

  IsoSpecWrapper::IsoSpecWrapper(double threshold, bool absolute) :
    threshold_(threshold),
    absolute_(absolute)
  {
  }

  std::vector<Peak1D> IsoSpecWrapper::run_(Iso& iso)
  {
    int tabSize = 1000;
    int hashSize = 1000;

    IsoThresholdGenerator generator(std::move(iso), threshold_, absolute_, tabSize, hashSize); 

    std::vector<Peak1D> distribution;

    while(generator.advanceToNextConfiguration())
        distribution.emplace_back(Peak1D(generator.mass(), generator.eprob()));

    return distribution;
  }

  std::vector<Peak1D> IsoSpecWrapper::run(const std::string& formula)
  {
    Iso iso(formula.c_str());
    return run_(iso);
  }

  std::vector<Peak1D> IsoSpecWrapper::run(const std::vector<int>& isotopeNr,
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
    const double** IM = new const double*[dimNumber];
    const double** IP = new const double*[dimNumber];
    for (int i=0; i<dimNumber; i++)
    {
      IM[i] = isotopeMasses[i].data();
      IP[i] = isotopeProbabilities[i].data();
    }

    Iso iso(dimNumber, isotopeNr.data(), atomCounts.data(), IM, IP);

    delete[] IM;
    delete[] IP;
    
    return run_(iso);
  }

#if 0
  void runLayered(const std::string& formula)
  {
    Iso* iso = new Iso(formula.c_str());

    double delta = -10;
    int tabSize = 1000;
    int hashSize = 1000;

    IsoLayeredGenerator* generator = new IsoLayeredGenerator(std::move(*iso), delta, tabSize, hashSize);
    Tabulator<IsoLayeredGenerator>* tabulator = new Tabulator<IsoLayeredGenerator>(generator, true, true, true, true); 

    delete generator;
    delete tabulator;
    delete iso;
  }

  void runOrdered(const std::string& formula)
  {
    Iso* iso = new Iso(formula.c_str());

    int tabSize = 1000;
    int hashSize = 1000;

    IsoOrderedGenerator* generator = new IsoOrderedGenerator(std::move(*iso), tabSize, hashSize);
    Tabulator<IsoOrderedGenerator>* tabulator = new Tabulator<IsoOrderedGenerator>(generator, true, true, true, true); 

    delete generator;
    delete tabulator;
    delete iso;
  }
#endif

}

