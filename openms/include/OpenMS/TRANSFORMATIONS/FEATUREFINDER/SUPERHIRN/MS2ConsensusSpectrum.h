// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_MS2CONSENSUSSPECTRUM_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_MS2CONSENSUSSPECTRUM_H

#include <OpenMS/CONCEPT/Types.h>

#include <string>
#include <map>

namespace OpenMS
{

  class OPENMS_DLLAPI MS2ConsensusSpectrum
  {


    ////////////////////////////////////////////////
    // declaration of the private members:

private:

    // stores the MS2 fragments:
    std::multimap<double, MS2Fragment> MS2FragmentPeaks;


protected:

    ////////////////////////////////////////////////
    // declaration of the public members:


    double startTR;
    double endTR;
    int z;
    int apexScan;
    int startScan;
    int endScan;


public:

    double precursorMZ;
    double TR;


    // mass to charge tolerance for MS2 trace level:
    static double MS2_MZ_TOLERANCE;

    // class destructor
    ~MS2ConsensusSpectrum();

    // class constructor
    MS2ConsensusSpectrum();
    MS2ConsensusSpectrum(MS2Fragment *);
    MS2ConsensusSpectrum(double iPrecursorMZ, double iTR, int iChrg, int iApexScan);

    // class copy constructor
    MS2ConsensusSpectrum(const MS2ConsensusSpectrum &);
    // class copy constructor
    MS2ConsensusSpectrum(const MS2ConsensusSpectrum *);


    //////////////////////////////////////////////////
    // overload operators:
    MS2ConsensusSpectrum & operator=(const MS2ConsensusSpectrum &);
    bool operator==(const MS2ConsensusSpectrum &);
    MS2ConsensusSpectrum & operator<=(const MS2ConsensusSpectrum &);
    MS2ConsensusSpectrum & operator>=(const MS2ConsensusSpectrum &);
    MS2ConsensusSpectrum & operator<(const MS2ConsensusSpectrum &);
    MS2ConsensusSpectrum & operator>(const MS2ConsensusSpectrum &);


    //////////////////////////////////////////////////////
    // compute the similarity of the elution shape of the
    // MS2 fragment to this MS2 consensus spectrum
    double getLCElutionPeakSimilarity(MS2Fragment *);


    // add a MS2 fragment:
    void addMS2Fragment(MS2Fragment *);

    // compute MS2 parameters
    void computeMS2SpectrumParameters();

    //*** PK removed, never used

    // process the stored fragments:
//  void processConsenusSpectraFragments();

    // remove outlier fragments based on their:
    // MS2Fragment::OutlierAttribute = ...
    // 1: retention time
    // 2: precursor mass
    // etc.
//  void removeOutlierFragments();

    // remove H2O loss region of the MS2 spectra
//  void removeWaterLossRegion( );


    // show MS2 spectrum info:
    void show_info();


    // find a corresponding MS2 fragment
    MS2Fragment * findMS2Fragment(double);



    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class

    // precursor mass:
    double getPrecursorMZ(){return precursorMZ; }

    // TR:
    double getTR(){return TR; }
    // start TR
    double getStartTR(){return startTR; }
    // end TR
    double getEndTR(){return endTR; }


    // set / get  the charge state of the precursor MZ:
    void setPrecursorChrg(int IN){ z = IN; }
    int getPrecursorChrg(){ return z; }
    // apex scan:
    int getApexScan(){return apexScan; }
    // start scan
    int getStartScan(){return startScan; }
    // end scan
    int getEndScan(){return endScan; }
    // get the number of consensus fragments:
    int getNbMS2Fragments(){return (int) MS2FragmentPeaks.size(); }

    // get the MS2 fragments list iterator:
    std::multimap<double, MS2Fragment>::iterator getMS2FragmentPeakStart(){return MS2FragmentPeaks.begin(); }
    std::multimap<double, MS2Fragment>::iterator getMS2FragmentPeakEnd(){return MS2FragmentPeaks.end(); }
    std::multimap<double, MS2Fragment> * getMS2FragmentMap(){return &MS2FragmentPeaks; }


  };

} // ns

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_MS2CONSENSUSSPECTRUM_H
