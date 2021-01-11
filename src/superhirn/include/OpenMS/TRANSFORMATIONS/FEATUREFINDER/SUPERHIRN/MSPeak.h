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
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//

#pragma once

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnConfig.h>

#include <string>
#include <vector>
#include <map>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>

namespace OpenMS
{

//class CentroidPeak; // changed to include

  class SUPERHIRN_DLLAPI MSPeak
  {

    ////////////////////////////////////////////////
    // declaration of the private members:

    double precursorMZ;
    double MZ;
    float INTENSITY;
    int SCAN;
    double TR;
    unsigned int CHRG;
    unsigned int NRISOTOPES;
    float SCORE;

    std::string extraMSPeakInfo;

    // child scan options:
    bool precursorMass;
    int childScan;

    double SignalToNoise;
    std::vector<CentroidPeak> ISOPEAKS;

private:

    ////////////////////////////////////////////////
    // declaration of the public members:

public:

    // class destructor
    ~MSPeak();

    // class constructor
    MSPeak(int, double, float);
    MSPeak();
    MSPeak(int, double, float, unsigned int, unsigned int, float, std::vector<CentroidPeak>);
    MSPeak(const MSPeak &);
    MSPeak(const MSPeak *);

    //////////////////////////////////////////////////
    // overload operators:
    MSPeak & operator=(const MSPeak &);
    MSPeak & operator<=(const MSPeak &);
    MSPeak & operator>=(const MSPeak &);
    MSPeak & operator<(const MSPeak &);
    MSPeak & operator>(const MSPeak &);

    // show content of peak:
    void show_info();

    // store the MS/MS scan number and activate this peak as precursor peak:
    void activateAsPrecursorPeak(int);

    //////////////////////////////////////////////////
    // check if the input mass matches one of the isotopic masses
    bool checkIsotopeBelongingAndAdjustMass(double, double);

    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class

    std::vector<CentroidPeak> & get_isotopic_peaks();

    std::vector<CentroidPeak>::iterator get_isotopic_peaks_start();

    std::vector<CentroidPeak>::iterator get_isotopic_peaks_end();

    void setExtraPeakInfo(std::string in);

    std::string getExtraPeakInfo();

    // precursor mass of the MS2 scan:
    void setPrecursorMZ(double in);

    double getPrecursorMZ();

    // precursor mass charge state:
    void setPrecursorCHRG(int in);

    int getPrecursorCHRG();

    // check if this peak has been determined as precursor:
    bool getPrecursorActivation();

    int get_Chrg();

    void set_Chrg(int z);

    int get_Scan();

    float get_intensity();

    double get_MZ();

    int get_scan_number();

    void set_retention_time(double in);

    double get_retention_time();

    unsigned int get_charge_state();

    unsigned int get_nr_isotopes();

    float get_score();

    double getSignalToNoise();

    void setSignalToNoise(double in);

  };

} // ns

