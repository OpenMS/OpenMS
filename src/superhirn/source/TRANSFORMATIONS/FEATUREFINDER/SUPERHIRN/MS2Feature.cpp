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
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Feature.h>

#include <cstdio>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Fragment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ClusteredMS2ConsensusSpectrum.h>

namespace OpenMS
{

////////////////////////////////////////////////
// constructor for the object MS2Feature:
  MS2Feature::MS2Feature(MS2Fragment * in) :
    ClusteredMS2ConsensusSpectrum(in)
  {
    ID = -1;
  }

////////////////////////////////////////////////
// constructor for the object MS2ConsensusSpectrum:
  MS2Feature::MS2Feature(double iPrecursorMZ, double iTR, int iChrg, int iApexScan) :
    ClusteredMS2ConsensusSpectrum(iPrecursorMZ, iTR, iChrg, iApexScan)
  {
    ID = -1;
  }

//////////////////////////////////////////////////
// class desctructor of MS2Feature
  MS2Feature::~MS2Feature()
  {
  }

//////////////////////////////////////////////////
// class copy constructor of MS2Feature
  MS2Feature::MS2Feature(const MS2Feature & tmp) :
    ClusteredMS2ConsensusSpectrum(tmp)
  {
    ID = tmp.ID;
  }

//////////////////////////////////////////////////
// class copy constructor of MS2Feature
  MS2Feature::MS2Feature(const MS2Feature * tmp) :
    ClusteredMS2ConsensusSpectrum(tmp)
  {
    ID = tmp->ID;
  }

/////////////////////////////////////////////
// show info
  void MS2Feature::show_info()
  {
    //printf("DELETED");
  }

}
