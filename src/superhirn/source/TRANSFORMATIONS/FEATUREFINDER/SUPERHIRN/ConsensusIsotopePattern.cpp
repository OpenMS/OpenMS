// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//


#include <map>
#include <vector>
#include <cmath>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ConsensusIsotopePattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnUtil.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>

namespace OpenMS
{

  using namespace std;

////////////////////////////////////////////////
// constructor for the object ConsensusIsotopePattern:
  ConsensusIsotopePattern::ConsensusIsotopePattern()
  {
  }

//////////////////////////////////////////////////
// class desctructor of ConsensusIsotopePattern
  ConsensusIsotopePattern::~ConsensusIsotopePattern()
  {

    isotopesTrace_.clear();
    mzIsotopesStDev_.clear();
    intensIsotopesStDev_.clear();
    rawIsotopes_.clear();

  }

//////////////////////////////////////////////////
// class copy constructor of ConsensusIsotopePattern
  ConsensusIsotopePattern::ConsensusIsotopePattern(const ConsensusIsotopePattern & tmp)
  {
    isotopesTrace_ = tmp.isotopesTrace_;
    mzIsotopesStDev_ = tmp.mzIsotopesStDev_;
    intensIsotopesStDev_ = tmp.intensIsotopesStDev_;
    rawIsotopes_ = tmp.rawIsotopes_;
  }

//////////////////////////////////////////////////
// class copy constructor of ConsensusIsotopePattern
  ConsensusIsotopePattern::ConsensusIsotopePattern(const ConsensusIsotopePattern * tmp)
  {
    isotopesTrace_ = tmp->isotopesTrace_;
    mzIsotopesStDev_ = tmp->mzIsotopesStDev_;
    intensIsotopesStDev_ = tmp->intensIsotopesStDev_;
    rawIsotopes_ = tmp->rawIsotopes_;
  }

//////////////////////////////////////////////////
// copy constructor:
  ConsensusIsotopePattern & ConsensusIsotopePattern::operator=(const ConsensusIsotopePattern & tmp)
  {
    isotopesTrace_ = tmp.isotopesTrace_;
    mzIsotopesStDev_ = tmp.mzIsotopesStDev_;
    intensIsotopesStDev_ = tmp.intensIsotopesStDev_;
    rawIsotopes_ = tmp.rawIsotopes_;
    return *this;
  }

/////////////////////////////////////////////////
// order an isotope trace in the correct cluster:
  void ConsensusIsotopePattern::addIsotopeTrace(double mz, double intens)
  {

    map<double, pair<vector<double>, vector<double> > >::iterator F = rawIsotopes_.lower_bound(mz);
    bool match = false;
    if (F != rawIsotopes_.end())
    {

      // compute the delta:
      if (SuperHirnUtil::compareMassValuesAtPPMLevel(mz, (*F).first, SuperHirnParameters::instance()->getMzTolPpm()))
      {
        (*F).second.first.push_back(mz);
        (*F).second.second.push_back(mz);
        match = true;
      }
      else if (F != rawIsotopes_.begin())
      {
        --F;
        if (SuperHirnUtil::compareMassValuesAtPPMLevel(mz, (*F).first, SuperHirnParameters::instance()->getMzTolPpm()))
        {
          (*F).second.first.push_back(mz);
          (*F).second.second.push_back(mz);
          match = true;
        }

      }

    }

    if (!match)
    {
      vector<double> mzTmp;
      mzTmp.push_back(mz);
      vector<double> intensTmp;
      intensTmp.push_back(intens);
      rawIsotopes_.insert(make_pair(mz, make_pair(mzTmp, intensTmp)));
    }

  }

/////////////////////////////////////////////////
// construc the consenus pattern:
  void ConsensusIsotopePattern::constructConsusPattern()
  {

    map<double, pair<vector<double>, vector<double> > >::iterator I = rawIsotopes_.begin();
    while (I != rawIsotopes_.end())
    {
      // condens a isotope peak trace:
      condensIsotopePattern(&(*I).second);
      ++I;
    }

  }

// copied from simple_math
  pair<double, double> simple_math_AVERAGE_and_STDEV(vector<double> * in)
  {

    double AVERAGE = 0;
    double STDEV = 0;

    if (in->empty())
    {
      return make_pair(AVERAGE, STDEV);
    }

    if (in->size() > 1)
    {
      vector<double>::iterator START = in->begin();
      while (START != in->end())
      {
        AVERAGE += (*START);
        ++START;
      }
      AVERAGE /= double(in->size());

      START = in->begin();
      while (START != in->end())
      {
        STDEV += ((AVERAGE - (*START)) * (AVERAGE - (*START)));
        ++START;
      }
      STDEV /= double(in->size());
      STDEV = sqrt(STDEV);
      return make_pair(AVERAGE, STDEV);
    }
    else
    {
      return make_pair((*in->begin()), 0.0);
    }
  }

//////////////////////////////////////////////////
// condens the pattern, make averge peaks from the traces:
  void ConsensusIsotopePattern::condensIsotopePattern(pair<vector<double>, vector<double> > * in)
  {

    // mz
    pair<double, double> mz = simple_math_AVERAGE_and_STDEV(&(in->first));
    // intens:
    pair<double, double> intens = simple_math_AVERAGE_and_STDEV(&(in->second));

    isotopesTrace_.insert(make_pair(mz.first, intens.first));
    mzIsotopesStDev_.push_back(mz.second);
    intensIsotopesStDev_.push_back(intens.second);

  }

}
