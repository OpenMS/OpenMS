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

#ifndef MS2_FEATURE_H
#define MS2_FEATURE_H

#include "MS2Fragment.h"
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ClusteredMS2ConsensusSpectrum.h>

namespace OpenMS
{

  class OPENMS_DLLAPI MS2Feature :
    public ClusteredMS2ConsensusSpectrum
  {


    using ClusteredMS2ConsensusSpectrum::operator=;

    ////////////////////////////////////////////////
    // declaration of the private members:

private:

    int ID;

    ////////////////////////////////////////////////
    // declaration of the public members:

public:


    // class destructor
    ~MS2Feature();
    // class constructor
    MS2Feature() {}
    MS2Feature(MS2Fragment *);
    MS2Feature(double iPrecursorMZ, double iTR, int iChrg, int iApexScan);

    // class copy constructor
    MS2Feature(const MS2Feature &);
    // class copy constructor
    MS2Feature(const MS2Feature *);


    //////////////////////////////////////////////////
    // overload operators:
    //MS2Feature& operator=(const MS2Feature&);
    bool operator==(const MS2Feature &);
    MS2Feature & operator<=(const MS2Feature &);
    MS2Feature & operator>=(const MS2Feature &);
    MS2Feature & operator<(const MS2Feature &);
    MS2Feature & operator>(const MS2Feature &);


    // show info
    void show_info();


    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class

    void setID(int in){ID = in; }
    int getID(){return ID; }


  };

} // ns

#endif
