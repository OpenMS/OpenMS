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
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Simon, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/MzIdentMLDOMHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/FORMAT/XTandemInfile.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/METADATA/PeptideHit.h>

#include <QtCore/QFile>
#include <QtCore/QDir>
#include <QtCore/QProcess>

#include <iostream>
#include <cmath>
#include <string>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_TopPerc TopPerc

  @brief TopPerc facilitates the input to, the call of and output integration of Percolator.
  Percolator (http://per-colator.com/) is a tool to apply semi-supervised learning for peptide
  identification from shotgun proteomics datasets.

  @experimental This tool is work in progress and usage and input requirements might change.

  <center>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ MSGF+\f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN ="center" ROWSPAN=1> @ref TOPP_IDFilter</td>
  <td VALIGN="middle" ALIGN ="center" ROWSPAN=1> @ref TOPP_IDMapper</td>
  </tr>
  </table>
  </center>

  <p>Percolator is search engine sensitive, i.e. it's input has to vary, depending on the search engine.</p>

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_TopPerc.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_TopPerc.html

  Percolator is written by Lukas Käll (http://per-colator.com/ Copyright Lukas Käll <lukas.kall@scilifelab.se>)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPPercolator :
  public TOPPBase
{
public:
  TOPPPercolator() :
    TOPPBase("TopPerc", "Facilitate input to Percolator and reintegrate.", false)
  {
  }

protected:
  void preparePIN(vector<PeptideIdentification>& peptide_ids, bool is_decoy, TextFile& txt, int minCharge, int maxCharge)
  {
    char out_sep = '\t';
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

            String spec_ref = it->getMetaValue("spectrum_id").toQString().toStdString();           //TODO consider other spectraIDFormats or keep index only in metavalue
            vector<String> scan_id;
            spec_ref.split("scan=", scan_id);

            int label = 1;
            String SpecId = "target_SII_";
            if (is_decoy)
            {
              SpecId = "decoy_SII_";
              label = -1;
            }

            SpecId += scan_id[1] + "_" + String(rank) + "_" + scan_id[1] + "_" + String(charge) + "_" + String(rank);

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
            int pepLen = hit->getSequence().toString().length();
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
            bool enzN = isEnz(peptide_without_modifications.at(0), peptide_without_modifications.at(2), getStringOption_("enzyme"));
            bool enzC = isEnz(peptide_without_modifications.at(peptide_without_modifications.size() - 3), peptide_without_modifications.at(peptide_without_modifications.size() - 1), getStringOption_("enzyme"));
            int enzInt = countEnzymatic(hit->getSequence().toUnmodifiedString(), getStringOption_("enzyme"));

            String peptide_with_modifications = aaBefore + string(".") + hit->getSequence().toString() + string(".") + aaAfter;
            String protein = hit->getPeptideEvidences().front().getProteinAccession();

            // One PeptideSpectrumHit with all its features
            String lis = SpecId + out_sep + String(label) + out_sep + scan_id[1] + out_sep + (String)rawScore + out_sep +
                         (String)denovoScore + out_sep + (String)scoreRatio + out_sep + (String)energy + out_sep + (String)ln_eval +
                         out_sep + (String)isotopeError + out_sep + (String)lnExplainedIonCurrentRatio + out_sep +
                         (String)lnNTermIonCurrentRatio + out_sep + (String)lnCTermIonCurrentRatio + out_sep + (String)lnMS2IonCurrent
                         + out_sep + (String)expMass + out_sep + (String)pepLen + out_sep + (String)dM + out_sep + (String)absdM + out_sep +
                         (String)meanErrorTop7 + out_sep + (String)sqMeanErrorTop7 + out_sep + (String)stdevErrorTop7 +
                         out_sep + String(ss.str()) + String(enzN) + out_sep + String(enzC) + out_sep + String(enzInt) + out_sep + 
                         peptide_with_modifications + out_sep + protein + out_sep;

            // peptide Spectrum Hit pushed to the output file
            txt.addLine(lis);
          }
        }
      }
    }
  }

  // Function taken from Enzyme.h from Percolator
  bool isEnz(const char& n, const char& c, string enz)
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

  // Function taken from Enzyme.h from Percolator
  size_t countEnzymatic(String peptide, string enz)
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
  double rescaleFragmentFeature(double featureValue, int NumMatchedMainIons)
  {
    // Rescale the fragment features to penalize features calculated by few ions
    int numMatchedIonLimit = 7;
    int numerator = (1 + numMatchedIonLimit) * (1 + numMatchedIonLimit);
    int denominator = (1 + (min)(NumMatchedMainIons, numMatchedIonLimit)) * (1 + (min)(NumMatchedMainIons, numMatchedIonLimit));
    return featureValue * ((double)numerator / denominator);
  }

  void registerOptionsAndFlags_()
  {
    registerInputFile_("percolator_executable", "<executable>", "", "Path to the percolator binary", true, false, ListUtils::create<String>("skipexists"));
    registerInputFile_("in_target", "<file>", "", "Input target file");
    registerInputFile_("in_decoy", "<file>", "", "Input decoy file");
    setValidFormats_("in_target", ListUtils::create<String>("mzid"));
    setValidFormats_("in_decoy", ListUtils::create<String>("mzid"));

    registerOutputFile_("out", "<file>", "", "Output file", true);
    registerStringOption_("enzyme", "<enzyme>", "trypsin", "Type of enzyme: no_enzyme,elastase,pepsin,proteinasek,thermolysin,chymotrypsin,lys-n,lys-c,arg-c,asp-n,glu-c,trypsin", false);

//    registerOutputFile_("r", "<file>", "out", "Output tab delimited results to a file instead of stdout", false);
    registerOutputFile_("X", "<file>", "", "path to file in xml-output format (pout). Default is: pout.tab", false);
    registerFlag_("e", "read xml-input format (pin) from standard input");
    registerFlag_("Z", "Include decoys (PSMs, peptides and/or proteins) in the xml-output. Only available if -X is used.");
    registerDoubleOption_("p", "<value>", 0.0, "Cpos, penalty for mistakes made on positive examples. Set by cross validation if not specified.", false);
    registerDoubleOption_("n", "<value>", 0.0, "Cneg, penalty for mistakes made on negative examples. Set by cross validation if not specified.", false);
    registerDoubleOption_("F", "<value>", 0.01, "False discovery rate threshold to define positive examples in training. Set by cross validation if 0. Default is 0.01.", false);
    registerDoubleOption_("t", "<value>", 0.01, "False discovery rate threshold for evaluating best cross validation result and the reported end result. Default is 0.01.", false);
    registerIntOption_("i", "<number>", 0, "Maximal number of iterations", false);
    registerFlag_("x", "Quicker execution by reduced internal cross-validation.");
    registerDoubleOption_("f", "<value>", 0.6, "Fraction of the negative data set to be used as train set when only providing one negative set, remaining examples will be used as test set. Set to 0.6 by default.", false);
    registerOutputFile_("J", "<file>", "", "Output the computed features to the given file in tab-delimited format. A file with the features with the given file name will be created", false);
    registerInputFile_("k", "<file>", "", "Input file given in the deprecated pin-xml format generated by e.g. sqt2pin with the -k option", false);
    registerOutputFile_("w", "<file>", "", "Output final weights to the given file", false);
    registerInputFile_("W", "<file>", "", "Read initial weights to the given file", false);
    registerStringOption_("V", "<featurename>", "", "The most informative feature given as the feature name, can be negated to indicate that a lower value is better.", false);
    registerIntOption_("v", "<level>", 2, "Set verbosity of output: 0=no processing info, 5=all, default is 2", false);
    registerFlag_("u", "Use unit normalization [0-1] instead of standard deviation normalization");
    registerFlag_("R", "Measure performance on test set each iteration");
    registerFlag_("O", "Override error check and do not fall back on default score vector in case of suspect score vector");
    registerIntOption_("S", "<value>", 1, "Setting seed of the random number generator. Default value is 1", false);
    registerFlag_("K", "Retention time features calculated as in Klammer et al.");
    registerFlag_("D", "Include description of correct features");
    registerOutputFile_("B", "<file>", "", "Output tab delimited results for decoys into a file", false);
    registerFlag_("U", "Do not remove redundant peptides, keep all PSMS and exclude peptide level probabilities.");
    registerFlag_("s", "skip validation of input file against xml schema");
    registerFlag_("A", "output protein level probabilities");
    registerDoubleOption_("a", "<value>", 0.0, "Probability with which a present protein emits an associated peptide (to be used jointly with the -A option). Set by grid search if not specified.", false);
    registerDoubleOption_("b", "<value>", 0.0, "Probability of the creation of a peptide from noise (to be used jointly with the -A option). Set by grid search if not specified", false);
    registerDoubleOption_("G", "<value>", 0.0, "Prior probability of that a protein is present in the sample ( to be used with the -A option). Set by grid search if not specified", false);
    registerFlag_("g", "treat ties as if it were one protein (Only valid if option -A is active).");
    registerFlag_("I", "use pi_0 value when calculating empirical q-values (no effect if option Q is activated) (Only valid if option -A is active).");
    registerFlag_("q", "output empirical q-values and p-values (from target-decoy analysis) (Only valid if option -A is active).");
    registerFlag_("N", "disactivates the grouping of proteins with similar connectivity, for example if proteins P1 and P2 have the same peptides matching both of them, P1 and P2 will not be grouped as one protein (Only valid if option -A is active).");
    registerFlag_("E", "Proteins graph will not be separated in sub-graphs (Only valid if option -A is active).");
    registerFlag_("C", "it does not prune peptides with a very low score (~0.0) which means that if a peptide with a very low score is matching two proteins, when we prune the peptide,it will be duplicated to generate two new protein groups (Only valid if option -A is active).");
    registerIntOption_("d", "<value>", 0, "Setting depth 0 or 1 or 2 from low depth to high depth(less computational time) of the grid search for the estimation Alpha,Beta and Gamma parameters for fido(Only valid if option -A is active). Default value is 0", false);
    registerStringOption_("P", "<value>", "random", "Define the text pattern to identify the decoy proteins and/or PSMs, set this up if the label that identifies the decoys in the database is not the default (by default : random) (Only valid if option -A  is active).", false);
    registerFlag_("T", "Reduce the tree of proteins (removing low scored proteins) in order to estimate alpha,beta and gamma faster.(Only valid if option -A is active).");
    registerFlag_("Y", "Use target decoy competition to compute peptide probabilities.(recommended when using -A).");
    registerFlag_("H", "Q-value threshold that will be used in the computation of the MSE and ROC AUC score in the grid search (recommended 0.05 for normal size datasets and 0.1 for big size datasets).(Only valid if option -A is active).");
    registerFlag_("fido-truncation", "Proteins with a very low score (< 0.001) will be truncated (assigned 0.0 probability).(Only valid if option -A is active)");
    registerFlag_("Q", "Uses protein group level inference, each cluster of proteins is either present or not, therefore when grouping proteins discard all possible combinations for each group.(Only valid if option -A is active and -N is inactive).");
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // general variables and data to perform TopPerc
    //-------------------------------------------------------------
    vector<PeptideIdentification> peptide_ids;
    vector<ProteinIdentification> protein_ids;
    vector<PeptideIdentification> peptide_ids_d;
    vector<ProteinIdentification> protein_ids_d;

    //-------------------------------------------------------------
    // parsing parameters and crashing if mandatory parameters are missing
    //-------------------------------------------------------------
    String inputfile_target_name = getStringOption_("in_target").toQString().toStdString();
    writeDebug_(String("Input file of target: ") + inputfile_target_name, 1);
    if (inputfile_target_name == "")
    {
      writeLog_("No target input file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    String inputfile_decoy_name = getStringOption_("in_decoy").toQString().toStdString();
    writeDebug_(String("Input file of decoy: ") + inputfile_decoy_name, 1);
    if (inputfile_decoy_name == "")
    {
      writeLog_("No decoy input file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    String percolator_executable(getStringOption_("percolator_executable"));
    writeDebug_(String("Path to the percolator: ") + percolator_executable, 1);
    if (percolator_executable == "") //TODO     TOPPBase::findExecutable
    {
      writeLog_("No path to percolator specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // get the file extension of the input files to start the correct converter
    vector<String> input_target_file;
    vector<String> input_decoy_file;
    inputfile_target_name.split('.', input_target_file);
    String data_target = input_target_file[input_target_file.size() - 1];
    inputfile_decoy_name.split('.', input_decoy_file);
    String data_decoy = input_decoy_file[input_decoy_file.size() - 1];

    TextFile txt;
    char out_sep = '\t';

    //get Information about database search
    String datab = "";

    // converter for MSGF+ & Mascot Files
    if (data_target == "mzid" && data_decoy == "mzid")
    {
      datab = "MSGF+";

      // TODO FOR FUTURE DEVELOPMENT: check out without explicit parameter setting if input file is target or decoy!!!
      // Both input files are read in
      MzIdentMLFile().load(inputfile_target_name, protein_ids, peptide_ids);
      MzIdentMLFile().load(inputfile_decoy_name, protein_ids_d, peptide_ids_d);
      LOG_INFO << "Using IDs from" << protein_ids.back().getSearchEngine() << endl;

      // Open File and check if the Identifier is MSGF+
      if (peptide_ids.front().getIdentifier() == "MS-GF+" && peptide_ids_d.front().getIdentifier() == "MS-GF+")
      {

        // Find out how many possible charges are available
        int maxCharge = 0;
        int minCharge = 10;
        for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
        {
          for (vector<PeptideHit>::const_iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
          {
            if (hit->getCharge() > maxCharge)
            {
              maxCharge = hit->getCharge();
            }
            if (hit->getCharge() < minCharge)
            {
              minCharge = hit->getCharge();
            }
          }
        }

        // Create String of the charges for the header of the tab file
        stringstream ss;
        ss << "Charge" << minCharge << ", ";
        for (int j = minCharge + 1; j < maxCharge + 1; j++)
        {
          ss << "Charge" << j << ",";
        }

        // Create header for the features
        string featureset = "SpecId, Label,ScanNr, RawScore, DeNovoScore,ScoreRatio, Energy,lnEValue,IsotopeError, lnExplainedIonCurrentRatio,lnNTermIonCurrentRatio,lnCTermIonCurrentRatio,lnMS2IonCurrent,Mass,PepLen,dM,absdM,MeanErrorTop7,sqMeanErrorTop7,StdevErrorTop7," + ss.str() + "enzN,enzC,enzInt,Peptide,Proteins";
        StringList txt_header0 = ListUtils::create<String>(featureset);
        txt.addLine(ListUtils::concatenate(txt_header0, out_sep));
        LOG_INFO << "consuming target file" << endl;
        // get all the features from the target file
        preparePIN(peptide_ids, false, txt, minCharge, maxCharge);
        LOG_INFO << "consuming decoy file" << endl;
        // get all the features from the decoy file
        preparePIN(peptide_ids_d, true, txt, minCharge, maxCharge);
      }
      else if (peptide_ids.front().getIdentifier() == "Mascot" && peptide_ids_d.front().getIdentifier() == "Mascot")
      {
        // TODO: Mascot Implementation
      }
    }
    // converter for XTandem-Files
    // TODO IN FUTURE DEVELOPMENT: IMPLEMENT MZID READER FOR XTANDEMFILES
    else if (data_target == "idXML" && data_decoy == "idXML")
    {
      datab = "XTANDEM";
      IdXMLFile file;
      IdXMLFile decoy_file;
      file.load(getStringOption_("in_target"), protein_ids, peptide_ids);
      decoy_file.load(getStringOption_("in_decoy"), protein_ids_d, peptide_ids_d);

      // Find out how many possible charges are available
      int maxCharge = 0;
      int minCharge = 10;

      for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
      {
        for (vector<PeptideHit>::const_iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          if (hit->getCharge() > maxCharge)
          {
            maxCharge = hit->getCharge();
          }
          if (hit->getCharge() < minCharge)
          {
            minCharge = hit->getCharge();
          }
        }
      }

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
          //TODO this must be spectrum_reference!!! parse spectrum number from there if necessary!
          String scannumber = String(it->getMetaValue("spectrum_id"));
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
          bool enzN = isEnz(peptide.at(0), peptide.at(2), getStringOption_("enzyme"));
          bool enzC = isEnz(peptide.at(peptide.size() - 3), peptide.at(peptide.size() - 1), getStringOption_("enzyme"));
          int enzInt = countEnzymatic(sequence, getStringOption_("enzyme"));
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
      for (vector<PeptideIdentification>::iterator it = peptide_ids_d.begin(); it != peptide_ids_d.end(); ++it)
      {
        if (it->isHigherScoreBetter())
        {
          String scannumber = String(it->getMetaValue("spectrum_id"));
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
          bool enzN = isEnz(peptide.at(0), peptide.at(2), getStringOption_("enzyme"));
          bool enzC = isEnz(peptide.at(peptide.size() - 3), peptide.at(peptide.size() - 1), getStringOption_("enzyme"));
          int enzInt = countEnzymatic(sequence, getStringOption_("enzyme"));
          String protein = it->getHits().front().getPeptideEvidences().front().getProteinAccession();

          // One PeptideSpectrumHit with all its features
          String lis = "_tandem_output_file_decoy_" + scannumber + "_" + String(charge) + "_1" + out_sep + String(label) + out_sep + scannumber + out_sep + String(hyperscore) + out_sep + String(deltascore) + out_sep + ss_ion_2.str() + out_sep
                       + String(mh) + out_sep + String(dm) + out_sep + String(absdM) + out_sep + String(length) + out_sep + ss.str() + out_sep + String(enzN) + out_sep + String(enzC) + out_sep + String(enzInt) + out_sep + peptide + out_sep + protein;

          // peptide Spectrum Hit pushed to the output file
          txt.addLine(lis);
        }
      }
    }
    else
    {
      LOG_INFO << "target and decoy files are not of the same type" << endl;
    }

    LOG_INFO << "Executing percolator" << endl;

    // create temp directory to store percolator in file pin.tab temporarily
    QDir qdir_temp(File::getTempDirectory().toQString());
    String temp_data_directory = File::getUniqueName();
    qdir_temp.mkdir(temp_data_directory.toQString());
    qdir_temp.cd(temp_data_directory.toQString());
    temp_data_directory  = File::getTempDirectory() + "/" + temp_data_directory;
    String in_file = temp_data_directory + "/" + File::getUniqueName() + ".tab";
    String out_file = temp_data_directory + "/" + File::getUniqueName() + ".tab";

    // File is stored in temp directory
    txt.store(in_file);

    QProcess process;
    QStringList arguments;

    // Check all set parameters and get them into arguments StringList
    arguments << "-r" << out_file.toQString();
    if (getFlag_("e")) arguments << "-e";
    if (getFlag_("Z")) arguments << "-Z";
    if (getDoubleOption_("p") != 0.0) arguments << "-p" << String(getDoubleOption_("p")).toQString();
    if (getDoubleOption_("n") != 0.0) arguments << "-n" << String(getDoubleOption_("n")).toQString();
    if (getDoubleOption_("F") != 0.01) arguments << "-F" << String(getDoubleOption_("F")).toQString();
    if (getDoubleOption_("t") != 0.01) arguments << "-t" << String(getDoubleOption_("t")).toQString();
    if (getIntOption_("i") != 0) arguments << "-i" << String(getIntOption_("i")).toQString();
    if (getFlag_("x")) arguments << "-x";
    if (getDoubleOption_("f") != 0.6) arguments << "-f" << String(getDoubleOption_("f")).toQString();
    if (getStringOption_("J") != "") arguments << "-J" << getStringOption_("J").toQString();
    if (getStringOption_("k") != "") arguments << "-k" << getStringOption_("k").toQString();
    if (getStringOption_("w") != "") arguments << "-w" << getStringOption_("w").toQString();
    if (getStringOption_("W") != "") arguments << "-W" << getStringOption_("W").toQString();
    if (getStringOption_("V") != "") arguments << "-V" << getStringOption_("V").toQString();
    if (getIntOption_("v") != 2) arguments << "-v" << String(getIntOption_("v")).toQString();
    if (getFlag_("u")) arguments << "-u";
    if (getFlag_("R")) arguments << "-R";
    if (getFlag_("O")) arguments << "-O";
    if (getIntOption_("S") != 1) arguments << "-S" << String(getDoubleOption_("S")).toQString();
    if (getFlag_("K")) arguments << "-K";
    if (getFlag_("D")) arguments << "-D";
    if (getStringOption_("B") != "") arguments << "-B" << getStringOption_("B").toQString();
    if (getFlag_("U")) arguments << "-U";
    if (getFlag_("s")) arguments << "-s";
    if (getFlag_("A")) arguments << "-A";
    if (getDoubleOption_("a") != 0.0) arguments << "-a" << String(getDoubleOption_("a")).toQString();
    if (getDoubleOption_("b") != 0.0) arguments << "-b" << String(getDoubleOption_("b")).toQString();
    if (getDoubleOption_("G") != 0.0) arguments << "-G" << String(getDoubleOption_("G")).toQString();
    if (getFlag_("g")) arguments << "-g";
    if (getFlag_("I")) arguments << "-I";
    if (getFlag_("q")) arguments << "-q";
    if (getFlag_("N")) arguments << "-N";
    if (getFlag_("E")) arguments << "-E";
    if (getFlag_("C")) arguments << "-C";
    if (getIntOption_("d") != 0) arguments << "-d" << String(getIntOption_("d")).toQString();
    if (getStringOption_("P") != "random") arguments << "-P" << getStringOption_("P").toQString();
    if (getFlag_("T")) arguments << "-T";
    if (getFlag_("Y")) arguments << "-Y";
    if (getFlag_("H")) arguments << "-H";
    if (getFlag_("fido-truncation")) arguments << "--fido-truncation";
    if (getFlag_("Q")) arguments << "-Q";
    arguments << "-U";
    arguments << in_file.toQString();

    // Percolator execution with the executable ant the arguments StringList
    process.execute(percolator_executable.toQString(), arguments); // does automatic escaping etc...

    // when percolator finished calculation, it stores the results -r option (with or without -U) or -m (which seems to be not working)
    CsvFile csv_file(out_file, '\t');

    map<String, vector<String> > pep_map;
    StringList row;

    for (UInt i = 1; i < csv_file.rowCount(); ++i)
    {
      csv_file.getRow(i, row);
      vector<String> row_values;
      // peptide
      row_values.push_back(row[4].chop(2).reverse().chop(2).reverse());
      // SVM-score
      row_values.push_back(row[1]);
      // Q-Value
      row_values.push_back(row[2]);
      // PEP
      row_values.push_back(row[3]);

      vector<String> substr;
      row[0].split('_', substr);
      pep_map[substr[2]] = row_values; // scannr. as written in preparePIN
    }


    // Add the percolator results to the peptide vector of the original input file
    for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
    {
      String sid = it->getMetaValue("spectrum_id");
      if (pep_map.find(sid) == pep_map.end())
      {
        vector<String> sr;
        sid.split('=', sr);
        sid = sr.back();
        if (pep_map.find(sid) == pep_map.end())
        {
          //no spectrum found - log?
          continue;
        }
      }

      it->setScoreType("q-value");
      it->setHigherScoreBetter(false);
      vector<PeptideHit> temp;
      swap(temp, it->getHits());
      for (vector<PeptideHit>::iterator hit = temp.begin(); hit != temp.end(); ++hit)
      {
        AASequence aat;
        aat.fromString(pep_map[sid][0]);
        if (hit->getSequence() == aat)
        {
          //get aa before/after/charge and metainfo
          hit->setMetaValue("MS:1001492", pep_map[sid][1].toDouble());        //svm score
          double qv = pep_map[sid][2].toDouble();        // q-value
          hit->setMetaValue("MS:1001491", qv);
          hit->setScore(qv);
          hit->setMetaValue("MS:1001493", pep_map[sid][3].toDouble());        //pep
          hit->setSequence(aat);
          it->insertHit(*hit);
        }
      }
      // TODO what with those not in percolator result file -> empty PeptideHit vector?
    }

    for (vector<ProteinIdentification>::iterator it = protein_ids.begin(); it != protein_ids.end(); ++it)
    {
      it->setSearchEngine("Percolator");
    }

    // Storing the PeptideHits with calculated q-value, pep and svm score
    MzIdentMLFile().store(getStringOption_("out").toQString().toStdString(), protein_ids, peptide_ids);

    LOG_INFO << "TopPerc finished successfully!" << endl;

    // As the percolator output file is not needed anymore, the temporary directory is going to be deleted
//    File::removeDirRecursively(temp_data_directory);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPPercolator tool;

  return tool.main(argc, argv);
}

/// @endcond
