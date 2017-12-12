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
// $Maintainer: Hannes Roest $
// $Authors: Witold Wolski, Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_DIAHELPER_H
#define OPENMS_ANALYSIS_OPENSWATH_DIAHELPER_H

#include <vector>
#include <OpenMS/CHEMISTRY/AASequence.h>

namespace OpenMS
{
  class TheoreticalSpectrumGenerator;
  namespace DIAHelpers
  {
    /**
      @brief Helper functions for the DIA scoring of OpenSWATH
    */
    ///@{
    /// compute the b and y series masses for a given AASequence
    OPENMS_DLLAPI void getBYSeries(AASequence& a,
                     std::vector<double>& bseries, std::vector<double>& yseries, TheoreticalSpectrumGenerator * g, UInt charge = 1u);

    /// for SWATH -- get the theoretical b and y series masses for a sequence
    OPENMS_DLLAPI void getTheorMasses(AASequence& a, std::vector<double>& masses, TheoreticalSpectrumGenerator * g, UInt charge = 1u);

    /// get averagine distribution given mass
    OPENMS_DLLAPI void getAveragineIsotopeDistribution(double product_mz,
                                         std::vector<std::pair<double, double> >& isotopesSpec, 
                                         double charge = 1.,
                                         int nr_isotopes = 4,
                                         double mannmass = 1.00048);

    /// simulate spectrum from AASequence
    OPENMS_DLLAPI void simulateSpectrumFromAASequence(AASequence& aa,
                                        std::vector<double>& firstIsotopeMasses, //[out]
                                        std::vector<std::pair<double, double> >& isotopeMasses, //[out]
                                        TheoreticalSpectrumGenerator * g,
                                        double charge = 1.);

    /// modify masses by charge
    OPENMS_DLLAPI void modifyMassesByCharge(const std::vector<std::pair<double, double> >& masses, //![in]
                              std::vector<std::pair<double, double> >& modmass, //!< [out]
                              double charge = 1.);

    /// add negative pre-isotope weights to spectrum
    OPENMS_DLLAPI void addPreisotopeWeights(const std::vector<double>& firstIsotopeMasses,
                              std::vector<std::pair<double, double> >& isotopeSpec, // output
                              UInt nrpeaks = 2, //nr of pre-isotope peaks
                              double preIsotopePeaksWeight = -0.5, // weight of pre-isotope peaks
                              double mannmass = 1.000482, //
                              double charge = 1.);

    /// given an experimental spectrum add isotope pattern.
    OPENMS_DLLAPI void addIsotopes2Spec(const std::vector<std::pair<double, double> >& spec,
                          std::vector<std::pair<double, double> >& isotopeMasses, //[out]
                          double charge = 1.);

    /// sorts vector of pairs by first
    OPENMS_DLLAPI void sortByFirst(std::vector<std::pair<double, double> >& tmp);
    /// extract first from vector of pairs
    OPENMS_DLLAPI void extractFirst(const std::vector<std::pair<double, double> >& peaks, std::vector<double>& mass);
    /// extract second from vector of pairs
    OPENMS_DLLAPI void extractSecond(const std::vector<std::pair<double, double> >& peaks, std::vector<double>& mass);
    
    ///}@
  }
}
#endif
