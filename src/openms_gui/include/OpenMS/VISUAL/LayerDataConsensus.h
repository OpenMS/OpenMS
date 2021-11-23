// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/LayerDataBase.h>

namespace OpenMS
{

  /**
  @brief Class that stores the data for one layer

  The data for a layer can be peak data, feature data (feature, consensus),
  chromatogram or peptide identification data. 

  For 2D and 3D data, the data is generally accessible through getPeakData()
  while features are accessible through getFeatureMap() and getConsensusMap().
  For 1D data, the current spectrum must be accessed through
  getCurrentSpectrum().

  Peak data is stored using a shared pointer to an MSExperiment data structure
  as well as a shared pointer to a OnDiscMSExperiment data structure. Note that
  the actual data may not be in memory as this is not efficient for large files
  and therefore may have to be retrieved from disk on-demand. 

  @note The spectrum for 1D viewing retrieved through getCurrentSpectrum() is a
  copy of the actual raw data and *different* from the one retrieved through
  getPeakData()[index]. Any changes to applied to getCurrentSpectrum() are
  non-persistent and will be gone the next time the cache is updated.
  Persistent changes can be applied to getPeakDataMuteable() and will be
  available on the next cache update.

  @note Layer is mainly used as a member variable of PlotCanvas which holds
  a vector of LayerData objects.

  @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI LayerDataConsensus : public LayerDataBase
  {
  public:
    /// Default constructor
    LayerDataConsensus() :
        LayerDataBase(LayerDataBase::DT_CONSENSUS){};
    /// no Copy-ctor (should not be needed)
    LayerDataConsensus(const LayerDataConsensus& ld) = delete;
    /// no assignment operator (should not be needed)
    LayerDataConsensus& operator=(const LayerDataConsensus& ld) = delete;
    /// move Ctor
    LayerDataConsensus(LayerDataConsensus&& ld) = default;
    /// move assignment
    LayerDataConsensus& operator=(LayerDataConsensus&& ld) = default;
  };

}// namespace OpenMS
