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
//  ( and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch)
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//


#pragma once

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnConfig.h>

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{

  class SUPERHIRN_DLLAPI MS2Fragment
  {


    ////////////////////////////////////////////////
    // declaration of the private members:

private:

    ////////////////////////////////////////////////
    // declaration of the public members:

    // AMRT tag
    double precursorMZ;
    int precursorCHRG;
    double TR;
    int scan;
    int z;

    double fragmentMZ;
    double intensityArea;

    // scan and TR ranges:
    int scanStart;
    int scanEnd;
    double trStart;
    double trEnd;


public:

    static int OutlierAttribute;

    // class destructor
    ~MS2Fragment();

    // constructor for the object MS2Fragment:
    MS2Fragment(double iPrecursorMZ, int iPrecursorCHRG, double iTR, int iScan, int iZ, double iFragmentMZ, double iIntensityArea,
                int iScanStart, int iScanEnd, double iTrStart, double iTrEnd);
    MS2Fragment(double iPrecursorMZ, int iPrecursorCHRG, double iTR, int iScan, int iZ, double iFragmentMZ, double iIntensityArea);


    // class constructor
    MS2Fragment();
    // class copy constructor
    MS2Fragment(const MS2Fragment &);
    // class copy constructor
    MS2Fragment(const MS2Fragment *);

    // show info of the MS2 fragment
    void show_info();

    // get the attribute of the fragment
    // according to which outliers are removed
    double getOutlierDetectionAttribute();


    //////////////////////////////////////////////////
    // overload operators:
    MS2Fragment & operator=(const MS2Fragment &);
    bool operator==(const MS2Fragment &);
    MS2Fragment & operator<=(const MS2Fragment &);
    MS2Fragment & operator>=(const MS2Fragment &);
    MS2Fragment & operator<(const MS2Fragment &);
    MS2Fragment & operator>(const MS2Fragment &);


    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class


    // get the averaged precursor mass:
    double getPrecursorMZ();
    void setPrecursorMZ(double iMZ);
    // get the averaged precursor charge:
    int getPrecursorCHRG();
    // retention time:
    double getTR();
    // start TR:
    double getStartTR();
    // end TR:
    double getEndTR();
    // get the Fragment MZ:
    double getFragmentMz();
    void setFragmentMz(double iMz);
    // get the charge state:
    int getCHRG();
    // get the apex scan:
    int getApexScan();
    // get the apex scan:
    int getStartScan();
    // get the apex scan:
    int getEndScan();

    // get the integrated peak area:
    double getFragmentPeakArea();
    void setFragmentPeakArea(double iIntens);

  };

} // ns

