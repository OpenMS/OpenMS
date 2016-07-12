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
    //TODO for all prepare* PSMId as written in PeptideIdentification::spectrum_reference
    //  and pre/post AA as - if begin/end of protein ([/] in PeptideEvidence) - see prepareMULTIpin
    //id <tab> label <tab> scannr <tab> feature1 <tab> ... <tab> featureN <tab> peptide <tab> proteinId1 <tab> .. <tab> proteinIdM

    void TopPerc::prepareCUSTOMpin(vector<PeptideIdentification>& peptide_ids, TextFile& txt, vector<String>& user_param_features, char out_sep)
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
          String scan_identifier = getScanIdentifier(it, peptide_ids.begin());
          Int scan_number = getScanNumber(scan_identifier);
          int label = 1;
          if (hit->metaValueExists("target_decoy") && String(hit->getMetaValue("target_decoy")).hasSubstring("decoy"))
          {
            label = -1;
          }

          StringList collected_feats;
          collected_feats.push_back(scan_identifier);
          collected_feats.push_back(String(label));
          collected_feats.push_back(String(scan_number));

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

    void TopPerc::prepareMSGFpin(vector<PeptideIdentification>& peptide_ids, string& enz, TextFile& txt, int min_charge, int max_charge, bool addMHC, char out_sep)
    {
      // Create String of the charges for the header of the tab file
      stringstream ss;
      ss << "Charge" << min_charge << ", ";
      for (int j = min_charge + 1; j < max_charge + 1; j++)
      {
        ss << "Charge" << j << ",";
      }

      // Create header for the features
      string featureset = "SpecId, Label,ScanNr, RawScore, DeNovoScore,ScoreRatio, Energy,lnEValue,IsotopeError, lnExplainedIonCurrentRatio,lnNTermIonCurrentRatio,lnCTermIonCurrentRatio,lnMS2IonCurrent,Mass,PepLen,dM,absdM,MeanErrorTop7,sqMeanErrorTop7,StdevErrorTop7," + ss.str() ;
      StringList txt_header = ListUtils::create<String>(featureset);
      if (addMHC)
      {
          txt_header.push_back("enzN");
          txt_header.push_back("enzC");
          txt_header.push_back("MHCLct");
          txt_header.push_back("Peptide");
          txt_header.push_back("Protein");
      }
      else
      {
          txt_header.push_back("enzN");
          txt_header.push_back("enzC");
          txt_header.push_back("enzInt");
          txt_header.push_back("Peptide");
          txt_header.push_back("Protein");
      }
      txt.addLine(ListUtils::concatenate(txt_header, out_sep));

      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        for (vector<PeptideHit>::const_iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          // Some Hits have no NumMatchedMainIons, and MeanError, etc. values. Have to ignore them!
          if (hit->metaValueExists("NumMatchedMainIons"))
          {
            // only take features from first ranked entries and only with meanerrortop7 != 0.0
            if (hit->getMetaValue("MeanErrorTop7").toString().toDouble() != 0.0)
            {
              int charge = hit->getCharge();

              String scan_identifier = getScanIdentifier(it, peptide_ids.begin());
              Int scan_number = getScanNumber(scan_identifier);
              int label = 1;
              if ((String(hit->getMetaValue("target_decoy"))).hasSubstring("decoy"))
              {
                label = -1;
              }

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
              int i = min_charge;
              while (i <= max_charge)
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
              String lis = scan_identifier + out_sep + String(label) + out_sep + String(scan_number) + out_sep + (String)rawScore + out_sep +
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

    void TopPerc::prepareXTANDEMpin(vector<PeptideIdentification>& peptide_ids, string& enz, TextFile& txt, int min_charge, int max_charge, char out_sep)
    {
      // Create String of the charges for the header of the tab file
      stringstream ss;
      ss << "Charge" << min_charge << ", ";
      for (int j = min_charge + 1; j < max_charge + 1; j++)
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
      StringList txt_header = ListUtils::create<String>(featureset);
      // Insert the header with the features names to the file
      txt.addLine(ListUtils::concatenate(txt_header, out_sep));

      LOG_INFO << "read in target file" << endl;
      // get all the features from the target file
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        if (it->isHigherScoreBetter())
        {
          String scan_identifier = getScanIdentifier(it, peptide_ids.begin());
          Int scan_number = getScanNumber(scan_identifier);
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
          int i = min_charge;
          while (i <= max_charge)
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
          String lis = scan_identifier + "_" + String(charge) +
            "_1" + out_sep + String(label) + out_sep + String(scan_number) + out_sep + String(hyperscore) +
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
          String scan_identifier = String(it->getMetaValue("spectrum_reference"));
          Int scan_number = getScanNumber(scan_identifier);
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
          int i = min_charge;
          while (i <= max_charge)
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
          String lis = scan_identifier + "_" + String(charge) + "_1" + out_sep + String(label) + out_sep + String(scan_number) + out_sep + String(hyperscore) + out_sep + String(deltascore) + out_sep + ss_ion_2.str() + out_sep
                       + String(mh) + out_sep + String(dm) + out_sep + String(absdM) + out_sep + String(length) + out_sep + ss.str() + out_sep + String(enzN) + out_sep + String(enzC) + out_sep + String(enzInt) + out_sep + peptide + out_sep + protein;

          // peptide Spectrum Hit pushed to the output file
          txt.addLine(lis);
        }
      }
    }

    void TopPerc::prepareCOMETpin(vector<PeptideIdentification>& peptide_ids, string& enz, TextFile& txt, int min_charge, int max_charge, char out_sep)
    {
      /** -with decoy comet search
id	label	ScanNr	lnrSp	deltLCn	deltCn	lnExpect	Xcorr	Sp	IonFrac	Mass	PepLen	Charge1	Charge2	Charge3	Charge4	Charge5	Charge6	enzN	enzC	enzInt	lnNumSP	dM	absdM	peptide	proteinId1
/home/.../150209_msms4_45_3_1	1	45	2.890372	0.066992	0.055908	2.212066	0.917335	94.189621	0.1875	1541.939549	13	0	0	1	0	0	0	0	0	2	9.352534	-0.000001	0.000001	H.FVIIIRKQTDLPV.I	XXX_sp|P30307|MPIP3_HUMAN
/home/.../150209_msms4_55_2_1	-1	55	2.70805	0.257442	0.142884	0.236914	0.884307	65.903358	0.3125	1087.697313	9	0	1	0	0	0	0	0	0	3	8.382976	0.000002	0.000002	F.TIRRKSLLT.S	DECOY_XXX_sp|Q5VYS8|TUT7_HUMAN
/home/.../150209_msms4_58_2_1	-1	58	2.079442	0.189541	0.038707	1.094758	0.910897	67.636864	0.2778	1086.695849	10	0	1	0	0	0	0	1	0	4	8.947416	-0.000003	0.000003	K.SKAKKPTKKA.K	DECOY_sp|P16403|H12_HUMAN
/home/.../150209_msms4_70_3_1	1	70	0.693147	0.199883	0.161329	0.813617	1.455887	249.314102	0.2708	1399.760349	13	0	0	1	0	0	0	0	1	1	10.987528	0.000009	0.000009	V.KFNGAHIPGSPFK.I	sp|Q14315|FLNC_HUMAN
      */

      // Create String of the charges for the header of the tab file
      stringstream ss;
      ss << "Charge" << min_charge << ", ";
      for (int j = min_charge+1; j <= max_charge; j++)
      {
        ss << "Charge" << j << ",";
      }

      String featureset = "id,label,ScanNr,lnrSp,deltLCn,deltCn,lnExpect,Xcorr,Sp,IonFrac,Mass,PepLen,"
              + ss.str()
              + "enzN,enzC,enzInt,lnNumSP,dM,absdM,peptide,proteinId1";
      StringList txt_header = ListUtils::create<String>(featureset);
      // Insert the header with the features names to the file
      txt.addLine(ListUtils::concatenate(txt_header, out_sep));

      // get all the feature values
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        double deltaLCn = 0;
        for (vector<PeptideHit>::iterator jt = it->getHits().begin(); jt != it->getHits().end(); ++jt)
        {
          deltaLCn += double(jt->getMetaValue("MS:1002253"));
        }
        it->sort();
        it->assignRanks();
        String scan_identifier = getScanIdentifier(it, peptide_ids.begin());
        Int scan_number = getScanNumber(scan_identifier);
        std::vector<PeptideHit> hits = it->getHits();
        for (vector<PeptideHit>::iterator jt = hits.begin(); jt != hits.end(); ++jt)
        {
          int charge = jt->getCharge();
          int label = 1;
          if (jt->metaValueExists("target_decoy") && String(jt->getMetaValue("target_decoy")).hasSubstring("decoy"))
          {
            label = -1;
          }
          //Xcorr
          String xcorr = String(jt->getMetaValue("MS:1002252"));
          //deltCn
          String deltaCn = String(jt->getMetaValue("MS:1002253"));
          //deltLCn deltaCn between first and last, i.e. sum in peptidehit
          //lnExpect
          String lnExpect = String(log(double(jt->getMetaValue("MS:1002257"))));
          //Sp
          String sp = String(jt->getMetaValue("MS:1002255"));
          //lnrSp log n rank Sp
          String lnrSp = String(log(double(jt->getMetaValue("MS:1002256"))));
          //IonFrac
          String ionfrac = double(jt->getMetaValue("MS:1002258"))/double(jt->getMetaValue("MS:1002259"));
          //Mass
          double mass = jt->getSequence().getMonoWeight(Residue::Full, charge)/charge;
          //PepLen
          int peplen = jt->getSequence().size(); //NB comet assigns peplen 0 to decoys?
          //Chargen
          StringList chargen;
          // write 1 for the correct charge, 0 for other charges
          for (int i = min_charge; i <= max_charge; ++i)
          {
             if (charge != i)
             {
               chargen.push_back("0");
             }
             else
             {
               chargen.push_back("1");
             }
          }
          //enzN
          bool enzN = isEnz(jt->getPeptideEvidences().front().getAABefore(), jt->getSequence().getPrefix(1).toString().c_str()[0], enz);
          //enzC
          bool enzC = isEnz(jt->getSequence().getSuffix(1).toString().c_str()[0], jt->getPeptideEvidences().front().getAAAfter(), enz);
          //enzInt
          int enzInt = countEnzymatic(jt->getSequence().toUnmodifiedString(), enz);
          //lnNumSP
          //this is practically not obtainable, as this seems to be the logn of the number of
          //internally matched decoy or target hits to that spectrum query depending on the current hit itself
          //is approximated by number of matched peptides
          String lnNumSP = String(log(double(jt->getMetaValue("num_matched_peptides"))));
          //dM
          double dm = it->getMZ() - mass;
          //absdM
          double absdm = abs(dm);
          //peptide
          String sequence = "";
          sequence += String(jt->getPeptideEvidences().front().getAABefore()); // just first peptide evidence
          sequence += jt->getSequence().toString();
          sequence += String(jt->getPeptideEvidences().front().getAAAfter()); //just first peptide evidence
          //proteinId1
          String pepevid = "";
          for (vector<PeptideEvidence>::const_iterator kt = jt->getPeptideEvidences().begin(); kt != jt->getPeptideEvidences().end(); ++kt)
          {
            pepevid += kt->getProteinAccession();
          }

          StringList row;
          row.push_back(scan_identifier);
          row.push_back(label);
          row.push_back(String(scan_number));
          row.push_back(lnrSp);
          row.push_back(deltaLCn);
          row.push_back(deltaCn);
          row.push_back(lnExpect);
          row.push_back(xcorr);
          row.push_back(sp);
          row.push_back(ionfrac);
          row.push_back(String(mass));
          row.push_back(String(peplen));
          row.push_back(ListUtils::concatenate(chargen, out_sep));
          row.push_back(String(enzN));
          row.push_back(String(enzC));
          row.push_back(String(enzInt));
          row.push_back(lnNumSP);
          row.push_back(String(dm));
          row.push_back(String(absdm));
          row.push_back(sequence);
          row.push_back(pepevid);

          txt.addLine(ListUtils::concatenate(row, out_sep));
        }
      }
    }

    void TopPerc::prepareMASCOTpin(vector<PeptideIdentification>& peptide_ids, string& enz, TextFile& txt, int min_charge, int max_charge, char out_sep)
    {
      /**
Features 1-9 Represent the Basic Feature Set and Features 10-18 Represent the Extended Feature Set As Used in Mascot Percolator

feature abbreviation	feature description
1. mass	Calculated monoisotopic mass of the identified peptide.
2. charge	Precursor ion charge
3. mScore	Mascot score
4. dScore	Mascot score minus Mascot score of next best nonisobaric peptide hit
5. deltaM	Calculated minus observed peptide mass (in Dalton and ppm).
6. absDeltaM	Absolute value of calculated minus observed peptide mass (in Dalton and ppm)
7. isoDeltaM	Calculated minus observed peptide mass, isotope error corrected (in Dalton and ppm)
8. uniquePeps	None (0), one (1), two or more (2) distinct peptide sequences match same protein
9. mc	Missed tryptic cleavages
10. totInt	Total ion intensity (log)
11. intMatchedTot	Total matched ion intensity (log)
12. relIntMatchedTot	Total matched ion intensity divided by total ion intensity
13. binom	Peptide Score as described in ref 28
14. fragMassError	Mean fragment mass error (in Dalton and ppm)
15. absFragMassError	Mean absolute fragment mass error (in Dalton and ppm)
16. fracIonsMatched	Fraction of calculated ions matched (per ion series)
17. seqCov	Sequence coverage of matched ions (per ion series)
18. intMatched	Matched ion intensity (per ion series)
      */
      // Create String of the charges for the header of the tab file
      stringstream ss;
      ss << "Charge" << min_charge << ", ";
      for (int j = min_charge+1; j <= max_charge; j++)
      {
        ss << "Charge" << j << ",";
      }

      String featureset = "id,label,ScanNr,mass,"
              + ss.str()
              + "mScore,dScore,deltaMass, absDeltaMass, uniqueToProt, enzN, enzC, enzInt, mod,sequence,protein";
      StringList txt_header = ListUtils::create<String>(featureset);
      // Insert the header with the features names to the file
      txt.addLine(ListUtils::concatenate(txt_header, out_sep));

      // get all the feature values
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        double deltaLCn = 0;
        for (vector<PeptideHit>::iterator jt = it->getHits().begin(); jt != it->getHits().end(); ++jt)
        {
          deltaLCn += double(jt->getMetaValue("MS:1002253"));
        }
        it->sort();
        it->assignRanks();
        String scan_identifier = getScanIdentifier(it, peptide_ids.begin());
        Int scan_number = getScanNumber(scan_identifier);
        it->sort();
        it->assignRanks();

        std::vector<PeptideHit> hits = it->getHits();
        assignDeltaScore(hits, "MS:1001171");
        for (vector<PeptideHit>::iterator jt = hits.begin(); jt != hits.end(); ++jt)
        {
          int label = 1;
          if (jt->metaValueExists("target_decoy") && String(jt->getMetaValue("target_decoy")).hasSubstring("decoy"))
          {
            label = -1;
          }
          int charge = jt->getCharge();
          double mass = jt->getSequence().getMonoWeight(Residue::Full, charge)/charge;
          //Chargen
          StringList chargen;
          // write 1 for the correct charge, 0 for other charges
          for (int i = min_charge; i <= max_charge; ++i)
          {
             if (charge != i)
             {
               chargen.push_back("0");
             }
             else
             {
               chargen.push_back("1");
             }
          }
          double mScore = double(jt->getMetaValue("MS:1001171"));
          double dScore = double(jt->getMetaValue("delta_score"));
          double dm = it->getMZ() - mass;
          double absdm = abs(dm);
          //no isoDeltaM - no isotope error info available from mascot adapter
          //no uniquePeps - no info from mascot substitute with sequence to protein uniqueness
          String uniquePeps = "0";
          if (String(jt->getMetaValue("protein_references")) == "unique")
          {
            uniquePeps = "1";
          }
          bool enzN = isEnz(jt->getPeptideEvidences().front().getAABefore(), jt->getSequence().getPrefix(1).toString().c_str()[0], enz);
          //enzC
          bool enzC = isEnz(jt->getSequence().getSuffix(1).toString().c_str()[0], jt->getPeptideEvidences().front().getAAAfter(), enz);
          //enzInt
          int enzInt = countEnzymatic(jt->getSequence().toUnmodifiedString(), enz);
          //no totInt info available from mascot adapter
          //no intMatchedTot info available from mascot adapter
          //no relIntMatchedTot info available from mascot adapter
          //no binom info available from mascot adapter
          //no fragMassError info available from mascot adapter
          //no absFragMassError	info available from mascot adapter
          //no fracIonsMatched info available from mascot adapter
          //no seqCov info available from mascot adapter
          //no intMatched info available from mascot adapter

          bool mod = jt->getSequence().isModified();
          String sequence = "";
          sequence += String(jt->getPeptideEvidences().front().getAABefore()); // just first peptide evidence
          sequence += jt->getSequence().toString();
          sequence += String(jt->getPeptideEvidences().front().getAAAfter()); //just first peptide evidence
          //proteinId1
          String pepevid = "";
          for (vector<PeptideEvidence>::const_iterator kt = jt->getPeptideEvidences().begin(); kt != jt->getPeptideEvidences().end(); ++kt)
          {
            pepevid += kt->getProteinAccession();
          }

          StringList row;
          row.push_back(scan_identifier);
          row.push_back(label);
          row.push_back(String(scan_number));
          row.push_back(String(mass));
          row.push_back(ListUtils::concatenate(chargen, out_sep));
          row.push_back(String(mScore));
          row.push_back(String(dScore));
          row.push_back(String(dm));
          row.push_back(String(absdm));
          row.push_back(String(uniquePeps));
          row.push_back(String(enzN));
          row.push_back(String(enzC));
          row.push_back(String(enzInt));
          row.push_back(String(mod));
          row.push_back(sequence);
          row.push_back(pepevid);

          txt.addLine(ListUtils::concatenate(row, out_sep));
        }
      }
    }

    // Function adapted from Enzyme.h in Percolator converter
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

    // Function adapted from Enzyme.h in Percolator converter
    Size TopPerc::countEnzymatic(String peptide, string& enz)
    {
      Size count = 0;
      for (Size ix = 1; ix < peptide.size(); ++ix)
      {
        if (isEnz(peptide[ix - 1], peptide[ix], enz))
        {
          ++count;
        }
      }
      return count;
    }

    // Function adapted from MsgfplusReader in Percolator converter
    double TopPerc::rescaleFragmentFeature(double featureValue, int NumMatchedMainIons)
    {
      // Rescale the fragment features to penalize features calculated by few ions
      int numMatchedIonLimit = 7;
      int numerator = (1 + numMatchedIonLimit) * (1 + numMatchedIonLimit);
      int denominator = (1 + (min)(NumMatchedMainIons, numMatchedIonLimit)) * (1 + (min)(NumMatchedMainIons, numMatchedIonLimit));
      return featureValue * ((double)numerator / denominator);
    }

    String TopPerc::getScanIdentifier(vector<PeptideIdentification>::iterator it, vector<PeptideIdentification>::iterator start)
    {
      String scan_identifier = it->getMetaValue("spectrum_reference");
      if (scan_identifier.empty())
      {
        scan_identifier = String(it->getMetaValue("spectrum_id"));
        if (scan_identifier.empty())
        {
          scan_identifier = String(it - start + 1);
          LOG_WARN << "no known spectrum identifiers, using index [1,n] - use at own risk." << endl;
        }
      }
      return scan_identifier.removeWhitespaces();
    }
    
    Int TopPerc::getScanNumber(String scan_identifier)
    {
      Size idx = 0;
      if ((idx = scan_identifier.find("index=")) != std::string::npos) 
      {
        scan_identifier = scan_identifier.substr(idx + 6);
      }
      else if ((idx = scan_identifier.find("scan=")) != std::string::npos) 
      {
        scan_identifier = scan_identifier.substr(idx + 5);
      }
      return scan_identifier.toInt();
    }

    void TopPerc::assignDeltaScore(vector<PeptideHit>& hits, String score_ref)
    {
      if (!hits.empty())
      {
        vector<PeptideHit>::iterator prev = hits.begin();
        double prev_score = double(prev->getMetaValue(score_ref));
        for (vector<PeptideHit>::iterator jt = hits.begin()+1; jt != hits.end(); ++jt)
        {
          double cur_score = double(jt->getMetaValue(score_ref));
          double value = prev_score - cur_score;
          prev->setMetaValue("delta_score",value);
          prev = jt;
        }
        (hits.end()-1)->setMetaValue("delta_score",0.0); //if last hit or only one hit
      }
    }

    void TopPerc::mergeMULTIids(vector<vector<ProteinIdentification> >& protein_ids_list, vector<vector<PeptideIdentification> >& peptide_ids_list, bool skip_checks)
    {
      //both input parameters must correspond
      if (peptide_ids_list.size() != protein_ids_list.size())
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Protein and Peptide Identification mismatch");
      }

      //search parameters of all runs must correspond (considering front() of each only)
      if (!skip_checks)
      {
        for (Size i=1; i < protein_ids_list.size(); ++i)
        {
          if( protein_ids_list[i-1].front().getSearchParameters().db != protein_ids_list[i].front().getSearchParameters().db )
          {
            throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, protein_ids_list[i-1].front().getSearchParameters().db+"!="+protein_ids_list[i].front().getSearchParameters().db);
          }
        }
      }

      LOG_DEBUG << "creating spectrum map" << endl;

      //setup map of merge characteristics per spectrum
      std::map<String,PeptideIdentification> unified;

      std::string common = "q-value_score, expect_score";
      StringList commonMetaValues = ListUtils::create<String>(common);
      for (std::vector<std::vector<PeptideIdentification> >::iterator pilit = peptide_ids_list.begin(); pilit != peptide_ids_list.end(); ++pilit)
      {
        String SE = protein_ids_list[distance(peptide_ids_list.begin(), pilit)].front().getSearchEngine();
        for (vector<PeptideIdentification>::iterator pit = pilit->begin(); pit != pilit->end(); ++pit)
        {
          PeptideIdentification ins = *pit;
          //prepare for merge
          for (vector<PeptideHit>::iterator hit = ins.getHits().begin(); hit != ins.getHits().end(); ++hit)
          {
            //move score from each hit to meta value
            hit->setMetaValue(SE + ":" + ins.getScoreType(), hit->getScore());
            //set score in each hit to #SE hits
            hit->setScore(1);
            //rename common meta values (to SE:commonmetavaluename)
            for (Size i = 0; i < commonMetaValues.size(); ++i)
            {
              if (hit->metaValueExists(commonMetaValues[i]))
              {
                DataValue val = hit->getMetaValue(commonMetaValues[i]);
                hit->setMetaValue(SE+":"+commonMetaValues[i],val);
                hit->removeMetaValue(commonMetaValues[i]);
              }
            }
          }
          ins.setScoreType("multiple");
          ins.setIdentifier("TopPerc_multiple_SE_input");
          String spectrum_reference = ins.getMetaValue("spectrum_reference");
          //merge in unified map
          if (unified.find(spectrum_reference) == unified.end())
          {
            unified[spectrum_reference] = ins;
          }
          else
          {
            //find corresponding hit
            for (vector<PeptideHit>::iterator hit = ins.getHits().begin(); hit != ins.getHits().end(); ++hit)
            {
              for (vector<PeptideHit>::iterator merger = unified[spectrum_reference].getHits().begin(); merger != unified[spectrum_reference].getHits().end(); ++merger)
              {
                if (hit->getSequence()==merger->getSequence())
                {
                  //care for peptide evidences!! set would be okay if checked for same search db in parameters,
//                  vector<PeptideEvidence> pev;
//                  pev.reserve(max(hit->getPeptideEvidences().size(),merger->getPeptideEvidences().size()));
//                  std::vector<ProteinHit>::iterator uni;
//                  std::sort(merger->getPeptideEvidences().begin(),merger->getPeptideEvidences().end(), TopPerc::lq_PeptideEvidence);
//                  std::sort(hit->getPeptideEvidences().begin(),hit->getPeptideEvidences().end(), TopPerc::lq_PeptideEvidence);
//                  uni = std::set_union(swop.front().getHits().begin(), swop.front().getHits().end(),
//                                       it->front().getHits().begin(),it->front().getHits().end(), pev.begin(),
//                                       TopPerc::lq_PeptideEvidence);
//                  pev.resize(uni-pev.begin());
//                  merger->setPeptideEvidences(pev);
                  //There is no mutable getPeptideEvidences() accessor in PeptideHit - above will not werk, but so long:
                  //Implying PeptideIndexer was applied (with the same search db each) will care for that all PeptideEvidences from two hits with equal AASequence are the same

                  //merge meta values
                  vector< String > keys;
                  hit->getKeys(keys);
                  for (vector<String>::const_iterator kt = keys.begin(); kt != keys.end(); ++kt)
                  {
                    if (!merger->metaValueExists(*kt))
                    {
                      merger->setMetaValue(*kt, hit->getMetaValue(*kt));
                    }
                  }
                  merger->setScore(merger->getScore() + hit->getScore());
                  break;
                }
              }
            }
          }
        }
      }
      LOG_DEBUG << "filled spectrum map" << endl;
      std::vector<PeptideIdentification> swip;
      swip.reserve(unified.size());
      LOG_DEBUG << "merging spectrum map" << endl;
      for (std::map<String,PeptideIdentification>::iterator it = unified.begin(); it != unified.end(); ++it)
      {
        swip.push_back(it->second);
      }
      peptide_ids_list.front().swap(swip);
      peptide_ids_list.resize(1);
      LOG_DEBUG << "Now containing " << peptide_ids_list.front().size() << " spectra identifications."<< endl;

      LOG_DEBUG << "merging search parameters" << endl;
      //care for search parameters!!

      std::vector<ProteinIdentification> swop;
      swop.push_back(ProteinIdentification());
      swop.back().setIdentifier("TopPerc_multiple_SE_input");
      swop.back().setDateTime(DateTime::currentDateTime());
      swop.back().setSearchParameters(protein_ids_list.front().front().getSearchParameters());
      for (std::vector<std::vector<ProteinIdentification> >::iterator it = protein_ids_list.begin(); it != protein_ids_list.end(); ++it)
      {
        std::vector<ProteinHit> v;
        v.resize(swop.front().getHits().size() + it->front().getHits().size());
        std::vector<ProteinHit>::iterator uni;
        std::sort(it->front().getHits().begin(),it->front().getHits().end(), TopPerc::lq_ProteinHit());
        LOG_DEBUG << "Sorted next part of the ProteinHits." << endl;
        LOG_DEBUG << "Melting with that many previous ProteinHits. " << swop.front().getHits().size() << endl;
        if (swop.front().getHits().empty())
        {
          v.swap(it->front().getHits());
        }
        else
        {
          uni = std::set_union(swop.front().getHits().begin(), swop.front().getHits().end(),
                             it->front().getHits().begin(),it->front().getHits().end(), v.begin(),
                             TopPerc::lq_ProteinHit());
          v.resize(uni-v.begin());
        }
        LOG_DEBUG << "Melting ProteinHits." << endl;

        swap(swop.front().getHits(),v);
        LOG_DEBUG << "Done with next ProteinHits." << endl;

        ProteinIdentification::SearchParameters sp = it->front().getSearchParameters();
        String SE = it->front().getSearchEngine();
        LOG_DEBUG << "Melting Parameters from " << SE << " into MetaInfo." << endl;
        {//insert into MetaInfo as SE:param
          swop.front().setMetaValue("SE:"+SE,it->front().getSearchEngineVersion());
          swop.front().setMetaValue(SE+":db",sp.db);
          swop.front().setMetaValue(SE+":db_version",sp.db_version);
          swop.front().setMetaValue(SE+":taxonomy",sp.taxonomy);
          swop.front().setMetaValue(SE+":charges",sp.charges);
          swop.front().setMetaValue(SE+":fixed_modifications",ListUtils::concatenate(sp.fixed_modifications, ","));
          swop.front().setMetaValue(SE+":variable_modifications",ListUtils::concatenate(sp.variable_modifications, ","));
          swop.front().setMetaValue(SE+":missed_cleavages",sp.missed_cleavages);
          swop.front().setMetaValue(SE+":fragment_mass_tolerance",sp.fragment_mass_tolerance);
          swop.front().setMetaValue(SE+":fragment_mass_tolerance_ppm",sp.fragment_mass_tolerance_ppm);
          swop.front().setMetaValue(SE+":precursor_tolerance",sp.precursor_tolerance);
          swop.front().setMetaValue(SE+":precursor_mass_tolerance_ppm",sp.precursor_mass_tolerance_ppm);
          swop.front().setMetaValue(SE+":digestion_enzyme",sp.digestion_enzyme.getName());
        }
        swop.front().setPrimaryMSRunPath(it->front().getPrimaryMSRunPath());
        swop.front().setSearchEngine("multiple");
        LOG_DEBUG << "Done with next Parameters." << endl;
      }
      protein_ids_list.front().swap(swop);
      protein_ids_list.resize(1);
      LOG_DEBUG << "All merging finished." << endl;

    }

    void TopPerc::prepareMULTIpin(vector<PeptideIdentification>& peptide_ids, ProteinIdentification& protein_id, string& enz, TextFile& txt, int min_charge, int max_charge, char out_sep)
    {
      //-------------------------------------------------------------
      // header
      //-------------------------------------------------------------
      // Create String of the charges for the header of the tab file
      stringstream ss;
      ss << "Charge" << min_charge << ", ";
      for (int j = min_charge+1; j <= max_charge; j++)
      {
        ss << "Charge" << j << ",";
      }

      StringList ses_used;
      StringList se_specifics;
      StringList keys;
      protein_id.getKeys(keys);

      if (ListUtils::contains(keys, "SE:MS-GF+"))
      {
        ses_used.push_back("MS-GF+");
        se_specifics.push_back("MS:1002049");  // rawscore
        se_specifics.push_back("MS:1002053");  // evalue
      }
      if (ListUtils::contains(keys, "SE:Mascot"))
      {
        ses_used.push_back("Mascot");
        se_specifics.push_back("Mascot_score");
        se_specifics.push_back("EValue");
      }
      if (ListUtils::contains(keys, "SE:Comet"))
      {
        ses_used.push_back("Comet");
        se_specifics.push_back("MS:1002252");  //xcorr
        se_specifics.push_back("MS:1002257");  //evalue
      }
      if (ListUtils::contains(keys, "SE:XTandem"))
      {
        ses_used.push_back("XTandem");
        se_specifics.push_back("XTandem_score");
        se_specifics.push_back("E-Value");
      }

      LOG_INFO << "Using " << ListUtils::concatenate(ses_used, ", ") << " as source for search engine specific features." << endl;

      String featureset = "id,label,ScanNr,"
              + ListUtils::concatenate(se_specifics, ",") + ","
              + ss.str()
              + "ionfrac,mass,enzN,enzC,enzInt,numHits,dM,absdM,PepLen,peptide,proteinId1";
      StringList txt_header = ListUtils::create<String>(featureset);
      // Insert the header with the features names to the file
      txt.addLine(ListUtils::concatenate(txt_header, out_sep));

      //-------------------------------------------------------------
      // values
      //-------------------------------------------------------------
      // get all the feature values
      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        it->sort();
        it->assignRanks();
        String scan_identifier = getScanIdentifier(it, peptide_ids.begin());
        Int scan_number = getScanNumber(scan_identifier);
        StringList idents;
        std::vector<PeptideHit> hits = it->getHits();
        for (vector<PeptideHit>::iterator jt = hits.begin(); jt != hits.end(); ++jt)
        {
          int charge = jt->getCharge();
          int label = 1;
          if (jt->metaValueExists("target_decoy") && String(jt->getMetaValue("target_decoy")).hasSubstring("decoy"))
          {
            label = -1;
          }

          StringList sesp;
          for (StringList::iterator s = se_specifics.begin(); s != se_specifics.end(); ++s)
          {
            if (jt->metaValueExists(*s))
              sesp.push_back(String(jt->getMetaValue(*s)));
            else
              sesp.push_back("-1");
          }

          StringList chargen;
          // write 1 for the correct charge, 0 for other charges
          for (int i = min_charge; i <= max_charge; ++i)
          {
             if (charge != i)
             {
               chargen.push_back("0");
             }
             else
             {
               chargen.push_back("1");
             }
          }

          //IonFrac
          String ionfrac = String(double(jt->getMetaValue("matched_intensity"))/double(jt->getMetaValue("sum_intensity")));  // also consider "matched_ion_number"/"peak_number"
          //Mass
          double mass = jt->getSequence().getMonoWeight(Residue::Full, charge)/charge;
          //enzN
          bool enzN = isEnz(jt->getPeptideEvidences().front().getAABefore(), jt->getSequence().getPrefix(1).toString().c_str()[0], enz);
          //enzC
          bool enzC = isEnz(jt->getSequence().getSuffix(1).toString().c_str()[0], jt->getPeptideEvidences().front().getAAAfter(), enz);
          //enzInt
          int enzInt = countEnzymatic(jt->getSequence().toUnmodifiedString(), enz);
          //numHits
          int numHits = jt->getScore();
          //dM
          double dm = it->getMZ() - mass;
          //absdM
          double absdm = abs(dm);
          //PepLen
          int peplen = jt->getSequence().size();
          //peptide
          String sequence = "";
          //replace flanking aa if [ or ] with -
          char pb = jt->getPeptideEvidences().front().getAABefore();
          sequence += pb=='['?"-.":String(pb)+"."; // just first peptide evidence
          sequence += jt->getSequence().toString();
          char pa = jt->getPeptideEvidences().front().getAAAfter();
          sequence += pa==']'?".-":"."+String(pa); // just first peptide evidence
          //proteinId1
          StringList pepevid;
          for (vector<PeptideEvidence>::const_iterator kt = jt->getPeptideEvidences().begin(); kt != jt->getPeptideEvidences().end(); ++kt)
          {
            pepevid.push_back(kt->getProteinAccession());
          }

          StringList row;
          row.push_back(scan_identifier);
          row.push_back(label);
          row.push_back(String(scan_number));
          row.push_back(ListUtils::concatenate(sesp, out_sep));
          row.push_back(ListUtils::concatenate(chargen, out_sep));
          row.push_back(ionfrac);
          row.push_back(String(mass));
          row.push_back(String(enzN));
          row.push_back(String(enzC));
          row.push_back(String(enzInt));
          row.push_back(String(numHits));
          row.push_back(String(dm));
          row.push_back(String(absdm));
          row.push_back(String(peplen));
          row.push_back(sequence);
          row.push_back(ListUtils::concatenate(pepevid, out_sep));

          txt.addLine(ListUtils::concatenate(row, out_sep));
        }
      }
    }

    void TopPerc::prepareCONCATpin(vector<vector<PeptideIdentification> >& peptide_id_list, vector<vector<ProteinIdentification> >& protein_id_list, string& enz, TextFile& txt, int min_charge, int max_charge, char out_sep)
    {
      //-------------------------------------------------------------
      // header
      //-------------------------------------------------------------
      // Create String of the charges for the header of the tab file
      stringstream ss;
      ss << "Charge" << min_charge << ", ";
      for (int j = min_charge+1; j <= max_charge; j++)
      {
        ss << "Charge" << j << ",";
      }

      StringList ses_used;
      for (vector<vector<ProteinIdentification> >::iterator it = protein_id_list.begin(); it != protein_id_list.end(); ++it)
      {
        ses_used.push_back(it->front().getSearchEngine());
      }

      LOG_INFO << "Using " << ListUtils::concatenate(ses_used, ", ") << " as source for search engine specific features." << endl;

      String featureset = "id,label,ScanNr,"
              + ListUtils::concatenate(ses_used, ",") + ","
              + ss.str()
              + "Evalue,ionfrac,mass,enzN,enzC,enzInt,numHits,dM,absdM,PepLen,peptide,proteinId1";
      StringList txt_header = ListUtils::create<String>(featureset);
      // Insert the header with the features names to the file
      txt.addLine(ListUtils::concatenate(txt_header, out_sep));

      //-------------------------------------------------------------
      // values
      //-------------------------------------------------------------
      // get all the feature values
      for (vector<vector<PeptideIdentification> >::iterator pit = peptide_id_list.begin(); pit != peptide_id_list.end(); ++pit)
      {
          Size i = std::distance(peptide_id_list.begin(),pit);
          String se = protein_id_list[i].front().getSearchEngine();
          for (vector<PeptideIdentification>::iterator it = pit->begin(); it != pit->end(); ++it)
          {
            it->sort();
            it->assignRanks();
            String scan_identifier = getScanIdentifier(it, pit->begin());
            Int scan_number = getScanNumber(scan_identifier);
            std::vector<PeptideHit> hits = it->getHits();
            for (vector<PeptideHit>::iterator jt = hits.begin(); jt != hits.end(); ++jt)
            {
              int charge = jt->getCharge();
              int label = 1;
              if (jt->metaValueExists("target_decoy") && String(jt->getMetaValue("target_decoy")).hasSubstring("decoy"))
              {
                label = -1;
              }

              StringList sesp;
              String ev;
              for (StringList::iterator s = ses_used.begin(); s != ses_used.end(); ++s)
              {
                if ((*s) == se)
                {
                  if (se == "MS-GF+")
                  {
                    sesp.push_back(jt->getMetaValue("MS:1002049"));  // rawscore
                    ev = jt->getMetaValue("MS:1002053");  // evalue
                  }
                  if (se == "Mascot")
                  {
                    sesp.push_back(jt->getMetaValue("Mascot_score"));
                    ev = jt->getMetaValue("EValue");
                  }
                  if (se == "Comet")
                  {
                    sesp.push_back(jt->getMetaValue("MS:1002252"));  //xcorr
                    ev = jt->getMetaValue("MS:1002257");  //evalue
                  }
                  if (se == "XTandem")
                  {
                    sesp.push_back(jt->getMetaValue("XTandem_score"));
                    ev = jt->getMetaValue("E-Value");
                  }
                }
                else
                  sesp.push_back("-1");
              }

              StringList chargen;
              // write 1 for the correct charge, 0 for other charges
              for (int i = min_charge; i <= max_charge; ++i)
              {
                 if (charge != i)
                 {
                   chargen.push_back("0");
                 }
                 else
                 {
                   chargen.push_back("1");
                 }
              }

              //IonFrac
              String ionfrac = String(double(jt->getMetaValue("matched_intensity"))/double(jt->getMetaValue("sum_intensity")));  // also consider "matched_ion_number"/"peak_number"
              //Mass
              double mass = jt->getSequence().getMonoWeight(Residue::Full, charge)/charge;
              //enzN
              bool enzN = isEnz(jt->getPeptideEvidences().front().getAABefore(), jt->getSequence().getPrefix(1).toString().c_str()[0], enz);
              //enzC
              bool enzC = isEnz(jt->getSequence().getSuffix(1).toString().c_str()[0], jt->getPeptideEvidences().front().getAAAfter(), enz);
              //enzInt
              int enzInt = countEnzymatic(jt->getSequence().toUnmodifiedString(), enz);
              //numHits
              int numHits = jt->getScore();
              //dM
              double dm = it->getMZ() - mass;
              //absdM
              double absdm = abs(dm);
              //PepLen
              int peplen = jt->getSequence().size();
              //peptide
              String sequence = "";
              //replace flanking aa if [ or ] with -
              char pb = jt->getPeptideEvidences().front().getAABefore();
              sequence += pb=='['?"-.":String(pb)+"."; // just first peptide evidence
              sequence += jt->getSequence().toString();
              char pa = jt->getPeptideEvidences().front().getAAAfter();
              sequence += pa==']'?".-":"."+String(pa); // just first peptide evidence
              //proteinId1
              StringList pepevid;
              for (vector<PeptideEvidence>::const_iterator kt = jt->getPeptideEvidences().begin(); kt != jt->getPeptideEvidences().end(); ++kt)
              {
                pepevid.push_back(kt->getProteinAccession());
              }

              StringList row;
              row.push_back(scan_identifier);
              row.push_back(label);
              row.push_back(String(scan_number));
              row.push_back(ListUtils::concatenate(sesp, out_sep));
              row.push_back(ListUtils::concatenate(chargen, out_sep));
              row.push_back(ev);
              row.push_back(ionfrac);
              row.push_back(String(mass));
              row.push_back(String(enzN));
              row.push_back(String(enzC));
              row.push_back(String(enzInt));
              row.push_back(String(numHits));
              row.push_back(String(dm));
              row.push_back(String(absdm));
              row.push_back(String(peplen));
              row.push_back(sequence);
              row.push_back(ListUtils::concatenate(pepevid, out_sep));

              txt.addLine(ListUtils::concatenate(row, out_sep));
          }
        }
      }
    }
    
    void TopPerc::readPoutAsMap(String pout_file, map<String, vector<PercolatorResult> >& pep_map)
    {
      CsvFile csv_file(pout_file, '\t');
      StringList row;

      for (Size i = 1; i < csv_file.rowCount(); ++i)
      {
        csv_file.getRow(i, row);
        PercolatorResult res(row);
        String spec_ref = res.PSMId;
        if (pep_map.find(spec_ref) == pep_map.end())
        {
          pep_map[spec_ref] = vector<TopPerc::PercolatorResult>();
        }
        pep_map[spec_ref].push_back(res);
      }
    }

}
