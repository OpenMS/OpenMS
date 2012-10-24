// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Florian Zeller $
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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_CONSENSISOTOPEPATTERN_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_CONSENSISOTOPEPATTERN_H

#include <OpenMS/CONCEPT/Types.h>

#include <map>
#include <vector>

namespace OpenMS
{

  class OPENMS_DLLAPI ConsensusIsotopePattern
  {

    ////////////////////////////////////////////////
    // declaration of the private members:

private:

    // stores the consensus pattern:
    std::map<double, double> isotopesTrace_;
    std::vector<double> mzIsotopesStDev_;
    std::vector<double> intensIsotopesStDev_;

    // stores the detected patterns by retention time
    std::map<double, std::pair<std::vector<double>, std::vector<double> > > rawIsotopes_;

    ////////////////////////////////////////////////
    // declaration of the public members:

public:

    // class destructor
    ~ConsensusIsotopePattern();

    // class constructor
    ConsensusIsotopePattern();
    // class copy constructor
    ConsensusIsotopePattern(const ConsensusIsotopePattern &);
    // class copy constructor
    ConsensusIsotopePattern(const ConsensusIsotopePattern *);

    //////////////////////////////////////////////////
    // overload operators:
    ConsensusIsotopePattern & operator=(const ConsensusIsotopePattern &);
    bool operator==(const ConsensusIsotopePattern &);
    ConsensusIsotopePattern & operator<=(const ConsensusIsotopePattern &);
    ConsensusIsotopePattern & operator>=(const ConsensusIsotopePattern &);
    ConsensusIsotopePattern & operator<(const ConsensusIsotopePattern &);
    ConsensusIsotopePattern & operator>(const ConsensusIsotopePattern &);

    // constructs the consensus pattern:
    void constructConsusPattern();
    // order an isotope trace in the correct cluster:
    void addIsotopeTrace(double, double);
    // condenses the pattern, make average peaks from the traces:
    void condensIsotopePattern(std::pair<std::vector<double>, std::vector<double> > *);

    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class

    std::map<double, double>::iterator getConsensIsotopeIteratorStart();
    std::map<double, double>::iterator getConsensIsotopeIteratorEnd();

  };

  inline std::map<double, double>::iterator ConsensusIsotopePattern::getConsensIsotopeIteratorStart()
  {
    return isotopesTrace_.begin();
  }

  inline std::map<double, double>::iterator ConsensusIsotopePattern::getConsensIsotopeIteratorEnd()
  {
    return isotopesTrace_.end();
  }

} // ns

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_CONSENSISOTOPEPATTERN_H
