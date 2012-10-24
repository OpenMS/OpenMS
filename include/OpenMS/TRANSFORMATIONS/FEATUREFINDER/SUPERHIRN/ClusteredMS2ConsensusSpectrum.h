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
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  ( and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch)
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_CLUSTEREDMS2CONSENSUSSPECTRUM_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_CLUSTEREDMS2CONSENSUSSPECTRUM_H

#include <vector>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2ConsensusSpectrum.h>

namespace OpenMS
{

  class OPENMS_DLLAPI ClusteredMS2ConsensusSpectrum :
    public MS2ConsensusSpectrum
  {

    ////////////////////////////////////////////////
    // declaration of the private members:

private:

    ////////////////////////////////////////////////
    // declaration of the public members:

    // stores the individual MS/MS spectra:
    std::vector<int> MS2Scans;

public:

    using MS2ConsensusSpectrum::operator=;

    // class destructor
    ~ClusteredMS2ConsensusSpectrum();

    // class constructor
    ClusteredMS2ConsensusSpectrum() {}
    ClusteredMS2ConsensusSpectrum(MS2Fragment *);
    ClusteredMS2ConsensusSpectrum(MS2ConsensusSpectrum *);
    ClusteredMS2ConsensusSpectrum(double, double, int, int);

    // class copy constructor
    ClusteredMS2ConsensusSpectrum(const ClusteredMS2ConsensusSpectrum &);
    // class copy constructor
    ClusteredMS2ConsensusSpectrum(const ClusteredMS2ConsensusSpectrum *);

    //////////////////////////////////////////////////
    // overload operators:
    //ClusteredMS2ConsensusSpectrum& operator=(const ClusteredMS2ConsensusSpectrum&);
    bool operator==(const ClusteredMS2ConsensusSpectrum &);
    ClusteredMS2ConsensusSpectrum & operator<=(const ClusteredMS2ConsensusSpectrum &);
    ClusteredMS2ConsensusSpectrum & operator>=(const ClusteredMS2ConsensusSpectrum &);
    ClusteredMS2ConsensusSpectrum & operator<(const ClusteredMS2ConsensusSpectrum &);
    ClusteredMS2ConsensusSpectrum & operator>(const ClusteredMS2ConsensusSpectrum &);

    //////////////////////////////////////////////////
    // trace the fragments across MS/MS scans runs:
    void constructClusteredConsenusSpectraFragments(MS2ConsensusSpectrum *);

    // add a MS2 Consensus Spectrum:
    void addMS2ConsensusSpectrum(MS2ConsensusSpectrum *);

    //////////////////////////////////////////////////
    // extracts fragments from a MS/MS spectra and inserts
    // them into the Clustered MS/MS spectrum:
    void extractFragmentsFromSpectra(MS2ConsensusSpectrum *);
    // merge a MS2 fragment into the target MS2 fragment:
    void mergeMS2Fragments(MS2Fragment *, MS2Fragment *);

    // plot all the consensus MS2 spectrum in one plot:
    void plotCombinedSpectra();

    // remove outlier fragments based on their:
    // MS2Fragment::OutlierAttribute = ...
    // 1: retention time
    // 2: precursor mass
    // etc.
    void removeOutlierFragments();

    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class

    int getNumberOfSpectraScan()
    {
      return (int) MS2Scans.size();
    }

    std::vector<int>::iterator getSpectraScanNumberStart()
    {
      return MS2Scans.begin();
    }

    std::vector<int>::iterator getSpectraScanNumberEnd()
    {
      return MS2Scans.end();
    }

  };

} // ns

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_CLUSTEREDMS2CONSENSUSSPECTRUM_H
