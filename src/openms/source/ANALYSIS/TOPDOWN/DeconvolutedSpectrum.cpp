// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------
//
// Created by Kyowon Jeong on 4/22/20.
//

#include "include/OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h"

namespace OpenMS
{
  DeconvolutedSpectrum::DeconvolutedSpectrum(MSSpectrum &s) :
      spec(&s)
  {
  }

  DeconvolutedSpectrum::~DeconvolutedSpectrum()
  {
    std::vector<PeakGroup>().swap(peakGroups);
  }

  bool DeconvolutedSpectrum::empty() const
  {
    return peakGroups.empty();
  }


  void DeconvolutedSpectrum::write(std::fstream &fs,
                                   FLASHDeconvHelperStructs::Parameter &param)//, fstream &fsm, fstream &fsp)
  {
    if (empty())
    {
      return;
    }

    for (auto &pg : peakGroups)
    {
      if (pg.peaks.empty())
      {
        continue;
      }
      double &m = pg.monoisotopicMass;
      double &am = pg.avgMass;
      double &intensity = pg.intensity;
      int minCharge = param.chargeRange + param.minCharge;
      int maxCharge = -1;
      for (auto &p : pg.peaks)
      {
        minCharge = minCharge < p.charge ? minCharge : p.charge;
        maxCharge = maxCharge > p.charge ? maxCharge : p.charge;
      }

      fs << pg.massIndex << "\t" << specIndex << "\t" << param.fileName << "\t" << spec->getNativeID() << "\t"
         << spec->getMSLevel() << "\t"
         << massCntr << "\t"
         << std::to_string(am) << "\t" << std::to_string(m) << "\t" << intensity << "\t"
         << (maxCharge - minCharge + 1) << "\t" << minCharge << "\t" << maxCharge << "\t"
         << std::to_string(spec->getRT())
         << "\t" << pg.peaks.size() << "\t";
      fs << std::fixed << std::setprecision(2);
      fs << pg.maxSNRcharge << "\t" << pg.maxSNR << "\t" << pg.maxSNRminMz << "\t" << pg.maxSNRmaxMz << "\t";
      fs << std::fixed << std::setprecision(-1);
      //if (n > 1 && pg.precursorScanNumber >= 0)
      //{
      // fs << pg.precursorSpecIndex << "\t" << pg.precursorMz << "\t" << pg.precursorCharge << "\t"
      //    << pg.precursorMonoMass << "\t" << pg.precursorIntensity << "\t" << pg.precursorSNR << "\t";
      //}

      if (param.writeDetail)
      {
        fs << std::fixed << std::setprecision(2);
        for (auto &p : pg.peaks)
        {
          fs << p.mz << ";";
        }
        fs << "\t";
        for (auto &p : pg.peaks)
        {
          fs << p.charge << ";";
        }
        fs << "\t";
        for (auto &p : pg.peaks)
        {
          fs << p.getUnchargedMass() << ";";
        }
        fs << "\t";
        for (auto &p : pg.peaks)
        {
          fs << p.isotopeIndex << ";";
        }
        fs << "\t";

        for (auto &p : pg.peaks)
        {
          auto tm = pg.monoisotopicMass + p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U;
          auto diff = (tm / p.charge + Constants::PROTON_MASS_U - p.mz) / p.mz;

          fs << 1e6 * diff << ";";
        }
        fs << "\t";

        fs << std::fixed << std::setprecision(1);
        for (auto &p : pg.peaks)
        {
          fs << p.intensity << ";";
        }
        fs << "\t";
      }

      fs << std::fixed << std::setprecision(3);
      fs << pg.isotopeCosineScore;

      if (spec->getMSLevel() == 1)
      {
        fs << "\t" << pg.chargeCosineScore;
      }

      fs << "\t" << pg.totalSNR << "\n" << std::setprecision(-1); // TODO
    }

  }


  void DeconvolutedSpectrum::writeHeader(std::fstream &fs, int &n, bool detail)
  {
    if (detail)
    {
      if (n == 1)
      {
        fs
            << "MassIndex\tSpecIndex\tFileName\tSpecID\tMSLevel\tMassCountInSpec\tAvgMass\tMonoisotopicMass\t"
               "AggregatedIntensity\tPeakChargeRange\tPeakMinCharge\tPeakMaxCharge\t"
               "RetentionTime\tPeakCount\tMaxSNRCharge\tMaxSNR\tMaxSNRMinMz\tMaxSNRMaxMz\t"
               //"PrecursorSpecIndex\tPrecursorMz\tPrecursorCharge\tPrecursorMonoMass\tPrecursorIntensity\t"
               "PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
               "PeakIntensities\tIsotopeCosineScore\tChargeIntensityCosineScore\tTotalSNR\tScore\n";
      }
      else
      {
        fs
            << "MassIndex\tSpecIndex\tFileName\tSpecID\tMSLevel\tMassCountInSpec\tAvgMass\tMonoisotopicMass\t"
               "AggregatedIntensity\tPeakChargeRange\tPeakMinCharge\tPeakMaxCharge\t"
               "RetentionTime\tPeakCount\tMaxSNRCharge\tMaxSNR\tMaxSNRMinMz\tMaxSNRMaxMz\t"
               "PrecursorSpecIndex\tPrecursorMz\tPrecursorCharge\tPrecursorMonoMass\tPrecursorIntensity\tPrecursorSNR\t"
               "PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
               "PeakIntensities\tIsotopeCosineScore\tTotalSNR\tScore\n";
      }
    }
    else
    {
      if (n == 1)
      {
        fs
            << "MassIndex\tSpecIndex\tFileName\tSpecID\tMSLevel\tMassCountInSpec\tAvgMass\tMonoisotopicMass\t"
               "AggregatedIntensity\tPeakChargeRange\tPeakMinCharge\tPeakMaxCharge\t"
               "RetentionTime\tPeakCount\tMaxSNRCharge\tMaxSNR\tMaxSNRMinMz\tMaxSNRMaxMz\t"
               //"PrecursorSpecIndex\tPrecursorMz\tPrecursorCharge\tPrecursorMonoMass\tPrecursorIntensity\t"
               //"PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakMzErrors\t"
               //"PeakIntensities\t"
               "IsotopeCosineScore\tChargeIntensityCosineScore\tTotalSNR\tScore\n";
      }
      else
      {
        fs
            << "MassIndex\tSpecIndex\tFileName\tSpecID\tMSLevel\tMassCountInSpec\tAvgMass\tMonoisotopicMass\t"
               "AggregatedIntensity\tPeakChargeRange\tPeakMinCharge\tPeakMaxCharge\t"
               "RetentionTime\tPeakCount\tMaxSNRCharge\tMaxSNR\tMaxSNRMinMz\tMaxSNRMaxMz\t"
               "PrecursorSpecIndex\tPrecursorMz\tPrecursorCharge\tPrecursorMonoMass\tPrecursorIntensity\tPrecursorSNR\t"
               //"PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakMzErrors\t"
               //"PeakIntensities\t"
               "IsotopeCosineScore\tTotalSNR\tScore\n";
      }

    }
    //pg.maxSNRcharge << "\t" << pg.maxSNR << "\t" << pg.maxSNRminMz << "\t" << pg.maxSNRmaxMz
    //MinScan	MaxScan	 RepScan	RepCharge	RepMz ApexScanNum    Envelope

  }


  void DeconvolutedSpectrum::writeTopFD(std::fstream &fs, int id)//, fstream &fsm, fstream &fsp)
  {
    if (spec->getMSLevel() == 1)
    {
      return;
    }
    fs << std::fixed << std::setprecision(2);
    fs << "BEGIN IONS\n"
       << "ID=" << id << "\n"
       << "SCANS=" << scanNumber << "\n"
       << "RETENTION_TIME=" << spec->getRT() << "\n";
    for (auto &a :  spec->getPrecursors()[0].getActivationMethods())
    {
      auto act = Precursor::NamesOfActivationMethodShort[a];
      fs << "ACTIVATION=" << Precursor::NamesOfActivationMethodShort[a] << "\n";
    }

    /*fs << "MS_ONE_ID=" << pg.precursorSpecIndex << "\n"
       << "MS_ONE_SCAN=" << pg.precursorScanNumber << "\n"
       << "PRECURSOR_MZ="
       << std::to_string(pg.precursorMz) << "\n"
       << "PRECURSOR_CHARGE=" << pg.precursorCharge << "\n"
       << "PRECURSOR_MASS=" << std::to_string(pg.precursorMonoMass) << "\n"
       << "PRECURSOR_INTENSITY=" << pg.precursorIntensity << "\n";*/
    fs << std::setprecision(-1);

    for (auto &pg : peakGroups)
    {
      fs << std::fixed << std::setprecision(2);
      fs << std::to_string(pg.monoisotopicMass) << "\t" << pg.intensity << "\t" << pg.maxSNRcharge
         //  << "\t" << log10(pg.precursorSNR+1e-10) << "\t" << log10(pg.precursorTotalSNR+1e-10)
         << "\t" <<
         log10(pg.maxSNR + 1e-10) << "\t" << log10(pg.totalSNR + 1e-10)
         << "\n";
      fs << std::setprecision(-1);
    }

    fs << "END IONS\n\n";
  }

}