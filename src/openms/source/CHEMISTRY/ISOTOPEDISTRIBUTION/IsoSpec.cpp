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

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsoSpec.h>

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <functional>
#include <numeric>
#include <iterator>

#include <OpenMS/CHEMISTRY/Element.h>
#include <include/OpenMS/CONCEPT/Constants.h>

#include "allocator.cpp"
#include "dirtyAllocator.cpp"
#include "isoSpec++.cpp"
#include "isoMath.cpp"
#include "marginalTrek++.cpp"
#include "operators.cpp"
#include "element_tables.cpp"
// #include "misc.cpp"
// #include "spectrum2.cpp" // contains unix-specific code
// #include "cwrapper.cpp"
#include "tabulator.cpp"

using namespace std;

/*
   void * setupIso(int             dimNumber,
   const int*      isotopeNumbers,
   const int*      atomCounts,
   const double*   isotopeMasses,
   const double*   isotopeProbabilities)
   */

namespace OpenMS
{

  IsoSpec::IsoSpec() :
    threshold_(0.01),
    absolute_(false)
  {
  }

  IsoSpec::IsoSpec(double threshold) :
    threshold_(threshold),
    absolute_(false)
  {
  }

  std::vector<double> IsoSpec::getMasses() {return masses_;}
  std::vector<double> IsoSpec::getProbabilities() {return probabilities_;}

  // Setup requires the following input:
  //    dimNumber = the number of elements (e.g. 3 for H, C, O)
  //    isotopeNumbers = a vector of how many isotopes each element has, e.g. [2, 2, 3])
  //    atomCounts = how many atoms of each we have [e.g. 12, 6, 6 for Glucose]
  //    isotopeMasses = array with a length of sum(isotopeNumbers) and the masses, e.g. [1.00782503227, 2.01410177819, 12, 13.0033548352, 15.9949146202, 16.9991317576, 17.9991596137]
  //    isotopeProbabilities = array with a length of sum(isotopeNumbers) and the probabilities, e.g. [0.999884, 0.0001157, 0.9892, 0.01078, etc ... ]
  //

#if 0
  void IsoSpec::run2(const std::string& formula)
  {
    Iso* iso = new Iso(formula.c_str());

    double delta = -10;
    int tabSize = 1000;
    int hashSize = 1000;

    IsoLayeredGenerator* generator = new IsoLayeredGenerator(std::move(*iso), delta, tabSize, hashSize);

    bool get_masses = true;
    bool get_probs = true;
    bool get_lprobs = true;
    bool get_confs = true;

    Tabulator<IsoLayeredGenerator>* tabulator = new Tabulator<IsoLayeredGenerator>(generator, get_masses, get_probs, get_lprobs, get_confs); 

    masses_.resize(tabulator->confs_no());
    masses_.assign(tabulator->masses(), tabulator->masses() + tabulator->confs_no());

    probabilities_.resize(tabulator->confs_no());
    probabilities_.assign(tabulator->probs(), tabulator->masses() + tabulator->confs_no());

    delete generator;
    delete tabulator;
  }

  void IsoSpec::run3(const std::string& formula)
  {
    Iso* iso = new Iso(formula.c_str());

    int tabSize = 1000;
    int hashSize = 1000;

    IsoOrderedGenerator* generator = new IsoOrderedGenerator(std::move(*iso), tabSize, hashSize);

    bool get_masses = true;
    bool get_probs = true;
    bool get_lprobs = true;
    bool get_confs = true;

    Tabulator<IsoOrderedGenerator>* tabulator = new Tabulator<IsoOrderedGenerator>(generator, get_masses, get_probs, get_lprobs, get_confs); 

    masses_.resize(tabulator->confs_no());
    masses_.assign(tabulator->masses(), tabulator->masses() + tabulator->confs_no());

    probabilities_.resize(tabulator->confs_no());
    probabilities_.assign(tabulator->probs(), tabulator->masses() + tabulator->confs_no());

    delete generator;
    delete tabulator;
  }

#endif

  void IsoSpec::run(const std::string& formula)
  {
    Iso* iso = new Iso(formula.c_str());

    int tabSize = 1000;
    int hashSize = 1000;

    IsoThresholdGenerator* generator = new IsoThresholdGenerator(std::move(*iso), threshold_, absolute_, tabSize, hashSize); 

    bool get_masses = true;
    bool get_probs = true;
    bool get_lprobs = true;
    bool get_confs = true;

    Tabulator<IsoThresholdGenerator>* tabulator = new Tabulator<IsoThresholdGenerator>(generator, get_masses, get_probs, get_lprobs, get_confs); 

    int size = tabulator->confs_no();

    masses_.clear();
    masses_.reserve(size);
    copy(&tabulator->masses()[0], &tabulator->masses()[size], back_inserter(masses_));

    probabilities_.clear();
    probabilities_.reserve(size);
    copy(&tabulator->probs()[0], &tabulator->probs()[size], back_inserter(probabilities_));

    delete generator;
    delete tabulator;
  }

  void run(const std::vector<int>& isotopeNr,
           const std::vector<int>& atomCounts,
           const std::vector<std::vector<double> >& isotopeMasses,
           const std::vector<std::vector<double> >& isotopeProbabilities)
  {
    int dimNumber = isotopeNr.size();
    // assert(isotopeNr.size() == atomCounts.size());
    // assert(isotopeMasses.size() == isotopeProbabilities.size());

    {

      const double** IM = new const double*[dimNumber];
      const double** IP = new const double*[dimNumber];
      for (int i=0; i<dimNumber; i++)
      {
          IM[i] = isotopeMasses[i].data();
          IP[i] = isotopeProbabilities[i].data();
      }

      Iso* iso = new Iso(dimNumber, isotopeNr.data(), atomCounts.data(), IM, IP);

      double threshold = 0.01;
      bool absolute = false;
      int tabSize = 1000;
      int hashSize = 1000;

      IsoThresholdGenerator* generator = new IsoThresholdGenerator(std::move(*iso), threshold, absolute, tabSize, hashSize); 

      bool get_masses = true;
      bool get_probs = true;
      bool get_lprobs = true;
      bool get_confs = true;

      Tabulator<IsoThresholdGenerator>* tabulator = new Tabulator<IsoThresholdGenerator>(generator, get_masses, get_probs, get_lprobs, get_confs); 

      const double* masses = tabulator->masses();
      const double* lprobs = tabulator->lprobs();
      const double* probs = tabulator->probs();
      const int* confs = tabulator->confs();

      delete generator;
      delete tabulator;
    }
  }


}

