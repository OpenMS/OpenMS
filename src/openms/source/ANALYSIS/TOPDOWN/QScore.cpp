//--------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include "OpenMS/ANALYSIS/TOPDOWN/QScore.h"
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <iomanip>

namespace OpenMS
{

  double QScore::getQScore(const PeakGroup *pg, const int abs_charge)
  {
    if (pg->empty())
    { // all zero
      return .0;
    }
    //const double th = 2;
    //const std::vector<double> weights_vh({1.3522, -1.0877, -16.4956, -2.036, -0.9439, 18.251});
    const std::vector<double> weights({ 1.492, -2.0041, -14.3891, -0.9853, 0.4568, 0.063, 14.4072});
    //const std::vector<double> weights({ -2.5334, -6.2718, -31.5409, 0.1305, -3.8576, 0.4948, 37.2419});

    //ChargeCos         1.492
    //ChargeSNR       -2.0041
    //Cos            -14.3891
    //SNR             -0.9853
    //ChargeScore      0.4568
    //AvgPPMerror       0.063
    //Intercept       14.4072

    // ChargeCos                    2.5334
    // ChargeSNR                    6.2718
    // Cos                         31.5409
    // SNR                         -0.1305
    // ChargeScore                  3.8576
    // AvgPPMerror                 -0.4948
    // Intercept                  -37.2419

    double score = weights.back();
    auto fv = toFeatureVector_(pg, abs_charge);

    for (Size i = 0; i < weights.size() - 1; i++)
    {
      score += fv[i] * weights[i];
    }
    double qscore = 1.0 / (1.0 + exp(score));
    /*if (qscore < th)
    {
      return qscore;
    }

    score = weights_h[weights_h.size() - 1];

    for (int i = 0; i < weights_h.size() - 1; i++)
    {
      score += fv[i] * weights_h[i];
    }
    qscore = 1.0 / (1.0 + exp(score));
*/
    return qscore;
  }

  std::vector<double> QScore::toFeatureVector_(const PeakGroup *pg, const int abs_charge)
  {
    std::vector<double> fvector(6);

    double a = pg->getChargeIsotopeCosine(abs_charge);
    double d = 1;
    int index = 0;
    fvector[index++] = (log2(a + d));
    a = pg->getChargeSNR(abs_charge);
    fvector[index++] = (log2(d + a / (1 + a)));
    a = pg->getIsotopeCosine();
    fvector[index++] = (log2(a + d));
    a = pg->getSNR();
    fvector[index++] = (log2(d + a / (1 + a)));
    a = pg->getChargeScore();
    fvector[index++] = (log2(a + d));
    a = pg->getAvgPPMError();
    fvector[index++] = (log2(abs(a*1e6) + d));
    return fvector;
  }

  void QScore::writeAttCsvFromTopPICHeader(std::fstream& f, bool write_detail)
  {
    f
        << "ACC,FirstResidue,LastResidue,ProID,RT,ScanNumber,PrecursorScanNumber,PrecursorMonoMass,PrecursorOriginalMonoMass,PrecursorAvgMass,PrecursorMz,PrecursorIntensity,"
           "MassIntensity,FeatureIntensity,PrecursorCharge,PrecursorMinCharge,"
           "PrecursorMaxCharge,PTM,PTMMass1,PTMMass2,PTMMass3,ChargeCos,ChargeSNR,Cos,SNR,ChargeScore,AvgPPMerror,Qscore,Evalue,Qvalue,";
    if (write_detail)
    {
      f << "PeakMZs,PeakIntensities,PeakMasses,PeakCharges,PeakIsotopeIndices,";
    }
    f << "Class\n";
  }

  void QScore::writeAttCsvFromTopPIC(const DeconvolvedSpectrum&  deconvolved_spectrum,
                           const FLASHDeconvHelperStructs::TopPicItem& top_id,
                           const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                           std::fstream& f,
                           bool write_detail)
  {
    int scan_number = deconvolved_spectrum.getScanNumber();
    double pmz = deconvolved_spectrum.getPrecursor().getMZ();
    auto pg = deconvolved_spectrum.getPrecursorPeakGroup();
    double pmass = //pg.getMonoMass();
        top_id.proteform_id < 0 ? pg.getMonoMass()
                                 : top_id.adj_precursor_mass;
    double precursor_intensity = deconvolved_spectrum.getPrecursor().getIntensity();
    int fr = top_id.first_residue;
    int lr = top_id.last_residue;
    String acc = top_id.protein_acc;
    int proID = top_id.proteform_id;
    double rt = deconvolved_spectrum.getOriginalSpectrum().getRT();
    double pscan = deconvolved_spectrum.getPrecursorScanNumber();
    double fintensity = top_id.intensity;
    int charge = deconvolved_spectrum.getPrecursorCharge();

    double e_value = top_id.e_value;
    double q_value = top_id.proteofrom_q_value;
    bool is_identified = top_id.proteform_id >= 0;
    auto ptm_mass = top_id.unexp_mod;

    auto avgpmass = avg.getAverageMassDelta(pmass) + pmass;
    if (pg.empty())
    {
      return;
    }
    else
    {
      auto fv = toFeatureVector_(&pg, charge);
      //if (pg.getChargeIsotopeCosine(charge) <= 0)
      //  return;
      double monomass = pmass;
      double mass = avgpmass;
      f << acc << "," << fr << "," << lr << "," << proID << "," << rt << "," << scan_number << "," << pscan << ","
        << monomass << "," << pg.getMonoMass() << "," << mass << "," << pmz << ","
        << precursor_intensity << ","
        << pg.getIntensity() << "," << fintensity << ","
        << charge << "," << std::get<0>(pg.getAbsChargeRange()) << "," << std::get<1>(pg.getAbsChargeRange()) << ","
        << (is_identified ? std::to_string(ptm_mass.size()) : "nan") << ",";
      for (Size k = 0; k < 3; k++)
      {
        if (k < ptm_mass.size())
        {
          f << ptm_mass[k] << ",";
        }
        else
        {
          f << "nan,";
        }
      }

      for (auto& item: fv)
      {
        f << item << ",";
      }

      f << pg.getQScore() << "," << e_value << "," << q_value << ",";
      if (write_detail)
      {
        f << std::fixed << std::setprecision(2);
        for (auto& p: pg)
        {
          f << p.mz << " ";
        }
        f << ";,";

        f << std::fixed << std::setprecision(1);
        for (auto& p: pg)
        {
          f << p.intensity << " ";
        }
        f << ";,";
        f << std::setprecision(-1);


        for (auto& p: pg)
        {
          f << p.getUnchargedMass() << " ";
        }
        f << ";,";

        for (auto& p: pg)
        {
          f << (p.is_positive ? p.abs_charge : -p.abs_charge) << " ";
        }
        f << ";,";

        for (auto& p: pg)
        {
          f << p.isotopeIndex << " ";
        }
        f << ";,";
        f << std::fixed << std::setprecision(-1);
      }


      f << (is_identified ? "T" : "F") << "\n";
    }
  }
  void QScore::writeAttCsvFromDecoyHeader(std::fstream& f)
  {
    f
      << "MSLevel,ChargeCos,ChargeSNR,Cos,SNR,ChargeScore,AvgPPMerror,Class\n";
    //ChargeCos         1.492
    //ChargeSNR       -2.0041
    //Cos            -14.3891
    //SNR             -0.9853
    //ChargeScore      0.4568
    //AvgPPMerror       0.063
  }
  void QScore::writeAttCsvFromDecoy(const DeconvolvedSpectrum& deconvolved_spectrum, std::fstream& f)
  {
    int ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
    for(auto& pg:deconvolved_spectrum)
    {
      String decoy =  pg.getDecoyIndex() == 0? "T" : "D";
      auto fv = toFeatureVector_(&pg, pg.getRepAbsCharge());
      f<< ms_level<<",";
      for (auto& item: fv)
      {
        f << item << ",";
      }
      f << decoy<< "\n";
    }
  }

}
