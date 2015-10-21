// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/ID/TopPerc.h>

using namespace std;

namespace OpenMS
{

    bool TopPerc::isEnz(const char& n, const char& c, string& enz)
    {
      if (enz == "trypsin")
      {
        return ((n == 'K' || n == 'R') && c != 'P') || n == '-' || c == '-';
      }
      else if (enz == "chymotrypsin")
      {
        return ((n == 'F' || n == 'W' || n == 'Y' || n == 'L') && c != 'P') || n == '-' || c == '-';
      }
      else if (enz == "thermolysin")
      {
        return ((c == 'A' || c == 'F' || c == 'I' || c == 'L' || c == 'M'
                || c == 'V' || (n == 'R' && c == 'G')) && n != 'D' && n != 'E') || n == '-' || c == '-';
      }
      else if (enz == "proteinasek")
      {
        return (n == 'A' || n == 'E' || n == 'F' || n == 'I' || n == 'L'
               || n == 'T' || n == 'V' || n == 'W' || n == 'Y') || n == '-' || c == '-';
      }
      else if (enz == "pepsin")
      {
        return ((c == 'F' || c == 'L' || c == 'W' || c == 'Y' || n == 'F'
                || n == 'L' || n == 'W' || n == 'Y') && n != 'R') || n == '-' || c == '-';
      }
      else if (enz == "elastase")
      {
        return ((n == 'L' || n == 'V' || n == 'A' || n == 'G') && c != 'P')
               || n == '-' || c == '-';
      }
      else if (enz == "lys-n")
      {
        return (c == 'K')
               || n == '-' || c == '-';
      }
      else if (enz == "lys-c")
      {
        return ((n == 'K') && c != 'P')
               || n == '-' || c == '-';
      }
      else if (enz == "arg-c")
      {
        return ((n == 'R') && c != 'P')
               || n == '-' || c == '-';
      }
      else if (enz == "asp-n")
      {
        return (c == 'D')
               || n == '-' || c == '-';
      }
      else if (enz == "glu-c")
      {
        return ((n == 'E') && (c != 'P'))
               || n == '-' || c == '-';
      }
      else
      {
        return true;
      }
    }

    void TopPerc::prepareCUSTOMpin(vector<PeptideIdentification>& peptide_ids, string& enz, TextFile& txt, vector<String>& user_param_features, char out_sep)
    {
      // Create header for the features
      string min_featureset = "SpecId, Label, ScanNr";
      StringList txt_header = ListUtils::create<String>(min_featureset);
      txt_header.insert(txt_header.end(), user_param_features.begin(), user_param_features.end() );
      txt.addLine(ListUtils::concatenate(txt_header, out_sep));

      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        for (vector<PeptideHit>::const_iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          String spec_ref = it->getMetaValue("spectrum_reference").toString();
          vector<String> scan_id;
          spec_ref.split("scan=", scan_id);
          String sid = scan_id.back();
          int label = 1;
          String SpecId = "target_SII_";
          if ((String(hit->getMetaValue("target_decoy"))).hasSubstring("decoy"))
          {
            SpecId = "decoy_SII_";
            label = -1;
          }

          SpecId += sid + "_" + String(hit->getCharge());

          StringList collected_feats;
          collected_feats.push_back(SpecId);
          collected_feats.push_back(label);
          collected_feats.push_back(sid);

          for (vector<String>::const_iterator feat = user_param_features.begin(); feat != user_param_features.end(); ++feat)
          {
          // Some Hits have no NumMatchedMainIons, and MeanError, etc. values. Have to ignore them!
            if (hit->metaValueExists(*feat))
            {
              collected_feats.push_back(hit->getMetaValue(*feat).toString());
            }
          }
          if (collected_feats.size() == user_param_features.size())
          { // only if all feats were present add
            txt.addLine(ListUtils::concatenate(collected_feats, out_sep));
          }
        }
      }
    }

    void TopPerc::prepareMSGFpin(vector<PeptideIdentification>& peptide_ids, string& enz, TextFile& txt, int minCharge, int maxCharge, bool addMHC, char out_sep)
    {
      // Create String of the charges for the header of the tab file
      stringstream ss;
      ss << "Charge" << minCharge << ", ";
      for (int j = minCharge + 1; j < maxCharge + 1; j++)
      {
        ss << "Charge" << j << ",";
      }

      // Create header for the features
      string featureset = "SpecId, Label,ScanNr, RawScore, DeNovoScore,ScoreRatio, Energy,lnEValue,IsotopeError, lnExplainedIonCurrentRatio,lnNTermIonCurrentRatio,lnCTermIonCurrentRatio,lnMS2IonCurrent,Mass,PepLen,dM,absdM,MeanErrorTop7,sqMeanErrorTop7,StdevErrorTop7," + ss.str() ;
      StringList txt_header0 = ListUtils::create<String>(featureset);
      if (addMHC)
      {
          txt_header0.push_back("enzN");
          txt_header0.push_back("enzC");
          txt_header0.push_back("MHCLct");
          txt_header0.push_back("Peptide");
          txt_header0.push_back("Protein");
      }
      else
      {
          txt_header0.push_back("enzN");
          txt_header0.push_back("enzC");
          txt_header0.push_back("enzInt");
          txt_header0.push_back("Peptide");
          txt_header0.push_back("Protein");
      }
      txt.addLine(ListUtils::concatenate(txt_header0, out_sep));

      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        for (vector<PeptideHit>::const_iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          // Some Hits have no NumMatchedMainIons, and MeanError, etc. values. Have to ignore them!
          if (hit->metaValueExists("NumMatchedMainIons"))
          {
            // only take features from first ranked entries and only with meanerrortop7 != 0.0
            if (hit->getRank() == 1 && hit->getMetaValue("MeanErrorTop7").toString().toDouble() != 0.0)
            {
              int rank = hit->getRank();
              int charge = hit->getCharge();

              String spec_ref = it->getMetaValue("spectrum_reference").toString();
              vector<String> scan_id;
              spec_ref.split("scan=", scan_id);
              String sid = scan_id.back();

              int label = 1;
              String SpecId = "target_SII_";
              if ((String(hit->getMetaValue("target_decoy"))).hasSubstring("decoy"))
              {
                SpecId = "decoy_SII_";
                label = -1;
              }

              SpecId += sid + "_" + String(rank) + "_" + sid + "_" + String(charge) + "_" + String(rank);

              double rawScore = hit->getMetaValue("MS:1002049").toString().toDouble();
              double denovoScore = hit->getMetaValue("MS:1002050").toString().toDouble();

              double scoreRatio;
              if (denovoScore > 0)
              {
                scoreRatio = (rawScore / denovoScore);
              }
              else
              {
                scoreRatio = rawScore * 10000;
              }

              double energy = denovoScore - rawScore;
              double ln_eval = -log(hit->getMetaValue("MS:1002053").toString().toDouble());
              int isotopeError = hit->getMetaValue("IsotopeError").toString().toInt();
              double lnExplainedIonCurrentRatio = log(hit->getMetaValue("ExplainedIonCurrentRatio").toString().toDouble() + 0.0001);           // @andsi: wtf?!
              double lnNTermIonCurrentRatio = log(hit->getMetaValue("NTermIonCurrentRatio").toString().toDouble() + 0.0001);           // @andsi: wtf?!
              double lnCTermIonCurrentRatio = log(hit->getMetaValue("CTermIonCurrentRatio").toString().toDouble() + 0.0001);           // @andsi: wtf?!
              double lnMS2IonCurrent = log(hit->getMetaValue("MS2IonCurrent").toString().toDouble());
              double expMass = it->getMZ();
              double calcMass = hit->getMetaValue("calcMZ");
              int pepLen = hit->getSequence().toUnmodifiedString().length();
              double dM = (expMass - (isotopeError * Constants::NEUTRON_MASS_U / charge) - calcMass) / expMass;
              double absdM = abs(dM);
              double meanErrorTop7 = hit->getMetaValue("MeanErrorTop7").toString().toDouble();
              int NumMatchedMainIons =  hit->getMetaValue("NumMatchedMainIons").toString().toInt();

              double stdevErrorTop7 = 0.0;
              if (hit->getMetaValue("StdevErrorTop7").toString() != "NaN")
              {
                stdevErrorTop7 = hit->getMetaValue("StdevErrorTop7").toString().toDouble();
                if (stdevErrorTop7 == 0.0)
                {
                  stdevErrorTop7 = meanErrorTop7;
                }
              }
              else
              {
                LOG_WARN << "Stdeverrortop7 is NaN" << endl;
              }

              meanErrorTop7 = rescaleFragmentFeature(meanErrorTop7, NumMatchedMainIons);
              double sqMeanErrorTop7 = rescaleFragmentFeature(meanErrorTop7 * meanErrorTop7, NumMatchedMainIons);
              stdevErrorTop7 = rescaleFragmentFeature(stdevErrorTop7, NumMatchedMainIons);

              // write 1 for the correct charge, 0 for other charges
              // i.e.: charge 3 for charges from 2-5: 0 1 0 0
              stringstream ss;
              int i = minCharge;
              while (i <= maxCharge)
              {
                if (charge != i)
                {
                  ss << "0" << out_sep;
                }
                if (charge == i)
                {
                  ss << "1" << out_sep;
                }
                i++;
              }
              char aaBefore = hit->getPeptideEvidences().front().getAABefore();
              char aaAfter = hit->getPeptideEvidences().front().getAAAfter();

              // sequence without modification: "ABC" instead of "ABC[UNIMOD:4]"
              String peptide_without_modifications = aaBefore + string(".") + hit->getSequence().toUnmodifiedString() + string(".") + aaAfter;

              // formula taken from percolator msgfplus-converter isEnz(n, c) for trypsin
              bool enzN = isEnz(peptide_without_modifications.at(0), peptide_without_modifications.at(2), enz);
              bool enzC = isEnz(peptide_without_modifications.at(peptide_without_modifications.size() - 3), peptide_without_modifications.at(peptide_without_modifications.size() - 1), enz);
              int enzInt = countEnzymatic(hit->getSequence().toUnmodifiedString(), enz);

              String peptide_with_modifications = aaBefore + string(".") + hit->getSequence().toString() + string(".") + aaAfter;
              String protein = hit->getPeptideEvidences().front().getProteinAccession();

              // One PeptideSpectrumHit with all its features
              String lis = SpecId + out_sep + String(label) + out_sep + scan_id[1] + out_sep + (String)rawScore + out_sep +
                           (String)denovoScore + out_sep + (String)scoreRatio + out_sep + (String)energy + out_sep + (String)ln_eval +
                           out_sep + (String)isotopeError + out_sep + (String)lnExplainedIonCurrentRatio + out_sep +
                           (String)lnNTermIonCurrentRatio + out_sep + (String)lnCTermIonCurrentRatio + out_sep + (String)lnMS2IonCurrent
                           + out_sep + (String)expMass + out_sep + (String)pepLen + out_sep + (String)dM + out_sep + (String)absdM + out_sep +
                           (String)meanErrorTop7 + out_sep + (String)sqMeanErrorTop7 + out_sep + (String)stdevErrorTop7 +
                           out_sep + String(ss.str());
              if (addMHC)
              {
                bool suf = false;
                static const string arr[] = {"A", "F", "I", "K", "M", "L", "R", "W", "V"};
                vector<string> mhcends (arr, arr + sizeof(arr) / sizeof(arr[0]) );
                for (std::vector<string>::iterator eit = mhcends.begin(); eit != mhcends.end(); ++eit)
                {
                  if (hit->getSequence().toUnmodifiedString().hasSuffix(string(*eit)))
                  {
                    suf = true;
                    break;
                  }
                }
                lis = lis + String(enzN) + out_sep + String(enzC) + out_sep
                        + String(suf) + out_sep + peptide_with_modifications + out_sep + protein + out_sep;
              }
              else
              {
                lis = lis + String(enzN) + out_sep + String(enzC) + out_sep
                        + String(enzInt) + out_sep + peptide_with_modifications + out_sep + protein + out_sep;
              }

              // peptide Spectrum Hit pushed to the output file
              txt.addLine(lis);
            }
          }
        }
      }
    }

    void TopPerc::prepareXTANDEMpin(vector<PeptideIdentification>& peptide_ids, string& enz, TextFile& txt, int minCharge, int maxCharge, char out_sep)
    {
      // Create String of the charges for the header of the tab file
      stringstream ss;
      ss << "Charge" << minCharge << ", ";
      for (int j = minCharge + 1; j < maxCharge + 1; j++)
      {

        ss << "Charge" << j << ",";
      }

      // Find out which ions are in XTandem-File and take only these as features
      stringstream ss_ion;
      if (peptide_ids.front().getHits().front().getMetaValue("a_score").toString() != "" &&
          peptide_ids.front().getHits().front().getMetaValue("a_ions").toString() != "")
      {
        ss_ion << "frac_ion_a" << ",";
      }
      if (peptide_ids.front().getHits().front().getMetaValue("b_score").toString() != "" &&
          peptide_ids.front().getHits().front().getMetaValue("b_ions").toString() != "")
      {
        ss_ion << "frac_ion_b" << ",";
      }
      if (peptide_ids.front().getHits().front().getMetaValue("c_score").toString() != "" &&
          peptide_ids.front().getHits().front().getMetaValue("c_ions").toString() != "")
      {
        ss_ion << "frac_ion_c" << ",";
      }
      if (peptide_ids.front().getHits().front().getMetaValue("x_score").toString() != "" &&
          peptide_ids.front().getHits().front().getMetaValue("x_ions").toString() != "")
      {
        ss_ion << "frac_ion_x" << ",";
      }
      if (peptide_ids.front().getHits().front().getMetaValue("y_score").toString() != "" &&
          peptide_ids.front().getHits().front().getMetaValue("y_ions").toString() != "")
      {
        ss_ion << "frac_ion_y" << ",";
      }
      if (peptide_ids.front().getHits().front().getMetaValue("z_score").toString() != "" &&
          peptide_ids.front().getHits().front().getMetaValue("z_ions").toString() != "")
      {
        ss_ion << "frac_ion_z" << ",";
      }

      // Create header for the features
      String featureset = "SpecId,Label,ScanNr,hyperscore,deltascore," + ss_ion.str() +
        ",Mass,dM,absdM,PepLen," + ss.str() + "enzN,enzC,enzInt,Peptide,Proteins";
      StringList txt_header0 = ListUtils::create<String>(featureset);
      // Insert the header with the features names to the file
      txt.addLine(ListUtils::concatenate(txt_header0, out_sep));

      LOG_INFO << "read in target file" << endl;
      // get all the features from the target file
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        if (it->isHigherScoreBetter())
        {
          String scannumber = String(it->getMetaValue("spectrum_reference"));
          int charge = it->getHits().front().getCharge();
          int label = 1;
          double hyperscore = it->getHits().front().getScore();
          // deltascore = hyperscore - nextscore
          double deltascore = hyperscore - it->getHits().front().getMetaValue("nextscore").toString().toDouble();
          String sequence = it->getHits().front().getSequence().toString();
          int length = sequence.length();

          // Find out correct ion types and get its Values
          stringstream ss_ion_2;

          if (it->getHits().front().getMetaValue("a_score").toString() != "" &&
              it->getHits().front().getMetaValue("a_ions").toString() != "")
          {
            ss_ion_2 << double(it->getHits().front().getMetaValue("a_ions")) / length << out_sep;
          }
          if (it->getHits().front().getMetaValue("b_score").toString() != "" &&
              it->getHits().front().getMetaValue("b_ions").toString() != "")
          {
            ss_ion_2 << double(it->getHits().front().getMetaValue("b_ions")) / length << out_sep;
          }
          if (it->getHits().front().getMetaValue("c_score").toString() != "" &&
              it->getHits().front().getMetaValue("c_ions").toString() != "")
          {
            ss_ion_2 << double(it->getHits().front().getMetaValue("c_ions")) / length << out_sep;
          }
          if (it->getHits().front().getMetaValue("x_score").toString() != "" &&
              it->getHits().front().getMetaValue("x_ions").toString() != "")
          {
            ss_ion_2 << double(it->getHits().front().getMetaValue("x_ions")) / length << out_sep;
          }
          if (it->getHits().front().getMetaValue("y_score").toString() != "" &&
              it->getHits().front().getMetaValue("y_ions").toString() != "")
          {
            ss_ion_2 << double(it->getHits().front().getMetaValue("y_ions")) / length << out_sep;
          }
          if (it->getHits().front().getMetaValue("z_score").toString() != "" &&
              it->getHits().front().getMetaValue("z_ions").toString() != "")
          {
            ss_ion_2 << double(it->getHits().front().getMetaValue("z_ions")) / length << out_sep;
          }
          double mass = it->getHits().front().getMetaValue("mass");
          double dm = it->getHits().front().getMetaValue("delta");
          double mh = mass + dm;
          double absdM = abs(dm);

          // write 1 for the correct charge, 0 for other charges
          // i.e.: charge 3 for charges from 2-5: 0 1 0 0
          stringstream ss;
          int i = minCharge;
          while (i <= maxCharge)
          {
            if (charge != i)
            {
              ss << "0" << out_sep;
            }
            if (charge == i)
            {
              ss << "1" << out_sep;
            }
            i++;
          }

          char aaBefore = it->getHits().front().getPeptideEvidences().front().getAABefore();
          char aaAfter = it->getHits().front().getPeptideEvidences().front().getAAAfter();

          String peptide = aaBefore + string(".") + sequence + string(".") + aaAfter;

          // formula taken from percolator converter isEnz(n, c) for trypsin
          bool enzN = isEnz(peptide.at(0), peptide.at(2), enz);
          bool enzC = isEnz(peptide.at(peptide.size() - 3), peptide.at(peptide.size() - 1), enz);
          int enzInt = countEnzymatic(sequence, enz);
          String protein = it->getHits().front().getPeptideEvidences().front().getProteinAccession();

          // One PeptideSpectrumHit with all its features
          String lis = "_tandem_output_file_target_" + scannumber + "_" + String(charge) +
            "_1" + out_sep + String(label) + out_sep + scannumber + out_sep + String(hyperscore) +
            out_sep + String(deltascore) + out_sep + ss_ion_2.str() + String(mh) + out_sep +
            String(dm) + out_sep + String(absdM) + out_sep + String(length) + out_sep + String(ss.str()) +
            String(enzN) + out_sep + String(enzC) + out_sep + String(enzInt) + out_sep + peptide + out_sep + protein;

          // peptide Spectrum Hit pushed to the output file
          txt.addLine(lis);
        }
      }

      LOG_INFO << "read in decoy file" << endl;
      // get all the features from the decoy file
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        if (it->isHigherScoreBetter())
        {
          String scannumber = String(it->getMetaValue("spectrum_reference"));
          int charge = it->getHits().front().getCharge();
          int label = -1;
          double hyperscore = it->getHits().front().getScore();
          // deltascore = hyperscore - nextscore
          double deltascore = hyperscore - it->getHits().front().getMetaValue("nextscore").toString().toDouble();
          String sequence = it->getHits().front().getSequence().toString();
          int length = sequence.length();

          // Find out correct ion types and get its Values
          stringstream ss_ion_2;

          if (it->getHits().front().getMetaValue("a_score").toString() != "" && it->getHits().front().getMetaValue("a_ions").toString() != "")
          {
            ss_ion_2 << double(it->getHits().front().getMetaValue("a_ions")) / length << out_sep;
          }
          if (it->getHits().front().getMetaValue("b_score").toString() != "" && it->getHits().front().getMetaValue("b_ions").toString() != "")
          {
            ss_ion_2 << double(it->getHits().front().getMetaValue("b_ions")) / length << out_sep;
          }
          if (it->getHits().front().getMetaValue("c_score").toString() != "" && it->getHits().front().getMetaValue("c_ions").toString() != "")
          {
            ss_ion_2 << double(it->getHits().front().getMetaValue("c_ions")) / length << out_sep;
          }
          if (it->getHits().front().getMetaValue("x_score").toString() != "" && it->getHits().front().getMetaValue("x_ions").toString() != "")
          {
            ss_ion_2 << double(it->getHits().front().getMetaValue("x_ions")) / length << out_sep;
          }
          if (it->getHits().front().getMetaValue("y_score").toString() != "" && it->getHits().front().getMetaValue("y_ions").toString() != "")
          {
            ss_ion_2 << double(it->getHits().front().getMetaValue("y_ions")) / length << out_sep;
          }
          if (it->getHits().front().getMetaValue("z_score").toString() != "" && it->getHits().front().getMetaValue("z_ions").toString() != "")
          {
            ss_ion_2 << double(it->getHits().front().getMetaValue("z_ions")) / length;
          }
          double mass = it->getHits().front().getMetaValue("mass");
          double dm = double(it->getHits().front().getMetaValue("delta"));
          double mh = mass + dm;
          double absdM = abs(dm);

          // write 1 for the correct charge, 0 for other charges
          // i.e: charge 3 for charges from 2-5: 0 1 0 0
          stringstream ss;
          int i = minCharge;
          while (i <= maxCharge)
          {
            if (charge != i)
            {
              ss << "0" << out_sep;
            }
            if (charge == i)
            {
              ss << "1" << out_sep;
            }
            i++;
          }

          char aaBefore = it->getHits().front().getPeptideEvidences().front().getAABefore();
          char aaAfter = it->getHits().front().getPeptideEvidences().front().getAAAfter();

          String peptide = aaBefore + string(".") + sequence + string(".") + aaAfter;

          // formula taken from percolator converter isEnz(n, c) for trypsin
          bool enzN = isEnz(peptide.at(0), peptide.at(2), enz);
          bool enzC = isEnz(peptide.at(peptide.size() - 3), peptide.at(peptide.size() - 1), enz);
          int enzInt = countEnzymatic(sequence, enz);
          String protein = it->getHits().front().getPeptideEvidences().front().getProteinAccession();

          // One PeptideSpectrumHit with all its features
          String lis = "_tandem_output_file_decoy_" + scannumber + "_" + String(charge) + "_1" + out_sep + String(label) + out_sep + scannumber + out_sep + String(hyperscore) + out_sep + String(deltascore) + out_sep + ss_ion_2.str() + out_sep
                       + String(mh) + out_sep + String(dm) + out_sep + String(absdM) + out_sep + String(length) + out_sep + ss.str() + out_sep + String(enzN) + out_sep + String(enzC) + out_sep + String(enzInt) + out_sep + peptide + out_sep + protein;

          // peptide Spectrum Hit pushed to the output file
          txt.addLine(lis);
        }
      }
    }

    // Function taken from Enzyme.h from Percolator
    size_t TopPerc::countEnzymatic(String peptide, string enz)
    {
      size_t count = 0;
      for (size_t ix = 1; ix < peptide.size(); ++ix)
      {
        if (isEnz(peptide[ix - 1], peptide[ix], enz))
        {
          ++count;
        }
      }
      return count;
    }

    // Function taken from the percolator converter MsgfplusReader
    double TopPerc::rescaleFragmentFeature(double featureValue, int NumMatchedMainIons)
    {
      // Rescale the fragment features to penalize features calculated by few ions
      int numMatchedIonLimit = 7;
      int numerator = (1 + numMatchedIonLimit) * (1 + numMatchedIonLimit);
      int denominator = (1 + (min)(NumMatchedMainIons, numMatchedIonLimit)) * (1 + (min)(NumMatchedMainIons, numMatchedIonLimit));
      return featureValue * ((double)numerator / denominator);
    }

}
