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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>

using namespace std;

namespace OpenMS
{
  GoodDiffFilter::GoodDiffFilter() :
    FilterFunctor()
  {
    setName(GoodDiffFilter::getProductName());
    //values from kinter sherman


    // TODO from CHEMISTRY!
    aamass_.insert(make_pair(57.02, 'G'));
    aamass_.insert(make_pair(71.04, 'A'));
    aamass_.insert(make_pair(87.03, 'S'));
    aamass_.insert(make_pair(97.05, 'P'));
    aamass_.insert(make_pair(99.07, 'V'));
    aamass_.insert(make_pair(101.05, 'T'));
    aamass_.insert(make_pair(103.01, 'C'));
    aamass_.insert(make_pair(113.08, 'L')); //and also I, but the chars are fillers, anyway
    aamass_.insert(make_pair(114.04, 'N'));
    aamass_.insert(make_pair(115.03, 'D'));
    aamass_.insert(make_pair(128.06, 'Q'));
    aamass_.insert(make_pair(128.09, 'K'));
    aamass_.insert(make_pair(129.04, 'E'));
    aamass_.insert(make_pair(131.04, 'M'));
    aamass_.insert(make_pair(137.06, 'H'));
    aamass_.insert(make_pair(147.07, 'F'));
    aamass_.insert(make_pair(156.10, 'R'));
    aamass_.insert(make_pair(163.06, 'Y'));
    aamass_.insert(make_pair(186.06, 'W'));

    //value from Bioinformatics, Bern 2004
    defaults_.setValue("tolerance", 0.37, "Tolerance value as defined by Bern et al.");
    defaultsToParam_();
  }

  GoodDiffFilter::GoodDiffFilter(const GoodDiffFilter & source) :
    FilterFunctor(source), aamass_(source.aamass_)
  {
  }

  GoodDiffFilter & GoodDiffFilter::operator=(const GoodDiffFilter & source)
  {
    if (this != &source)
    {
      FilterFunctor::operator=(source);
      aamass_ = source.aamass_;
    }
    return *this;
  }

  GoodDiffFilter::~GoodDiffFilter()
  {
  }

}
