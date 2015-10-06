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
// $Author: Mathias Walzer, Sven Nahnsen $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/QcMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/SYSTEM/File.h>

#include <QFileInfo>
//~ #include <boost/regex.hpp>

#include <vector>
#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_QCCalculator QCCalculator

    @brief Calculates basic quality parameters from MS experiments and compiles data for subsequent QC into a qcML file.

    <CENTER>
      <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ QCCalculator \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_QCMerger </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_XTandemAdapter </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_QCExporter </td>
        </tr>
      </table>
    </CENTER>

    The calculated quality parameters or data compiled as attachments for easy plotting input include file origin, spectra distribution, aquisition details, ion current stability ( & TIC ), id accuracy statistics and feature statistics.
    The MS experiments base name is used as name to the qcML element that is comprising all quality parameter values for the given run (including the given downstream analysis data).
    
    - @p id produces quality parameter values for the identification file;
    - @p feature produces quality parameter values for the feature file;
    - @p consensus produces quality parameter values for the consensus file;
    some quality parameter calculation are only available if both feature and ids are given.
    - @p remove_duplicate_features only needed when you work with a set of merged features. Then considers duplicate features only once.

    Output is in qcML format (see parameter @p out) which can be viewed directly in a modern browser (chromium, firefox, safari).

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_QCCalculator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_QCCalculator.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

class TOPPQCCalculator :
  public TOPPBase
{
public:
  TOPPQCCalculator() :
    TOPPBase("QCCalculator", "Calculates basic quality parameters from MS experiments and subsequent analysis data as identification or feature detection.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "raw data input file (this is relevant if you want to look at MS1, MS2 and precursor peak information)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Your qcML file.");
    setValidFormats_("out", ListUtils::create<String>("qcML"));
    registerInputFile_("id", "<file>", "", "Input idXML file containing the identifications. Your identifications will be exported in an easy-to-read format", false);
    setValidFormats_("id", ListUtils::create<String>("idXML"));
    registerInputFile_("feature", "<file>", "", "feature input file (this is relevant for most QC issues)", false);
    setValidFormats_("feature", ListUtils::create<String>("featureXML"));
    registerInputFile_("consensus", "<file>", "", "consensus input file (this is only used for charge state deconvoluted output. Use the consensusXML output form the DeCharger)", false);
    setValidFormats_("consensus", ListUtils::create<String>("consensusXML"));
    registerFlag_("remove_duplicate_features", "This flag should be set, if you work with a set of merged features.");
    //~ registerFlag_("MS1", "This flag should be set, if you want to work with MS1 stats.");
    //~ registerFlag_("MS2", "This flag should be set, if you want to work with MS2 stats.");
  }

  double getMassDifference(double theo_mz, double exp_mz, bool use_ppm)
  {
    double error(exp_mz - theo_mz);
    if (use_ppm)
    {
      error = error / (theo_mz * (double)1e-6);
      //~ error = (1-exp_mz/theo_mz) * (double)1e6;
    }
    return error;
  }

  ExitCodes main_(int, const char**)
  {
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    ProteinHit temp_protein_hit;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_id               = getStringOption_("id");
    String inputfile_feature       = getStringOption_("feature");
    String inputfile_consensus  = getStringOption_("consensus");
    String inputfile_raw            = getStringOption_("in");
    String outputfile_name       = getStringOption_("out");

    //~ bool Ms1(getFlag_("MS1"));
    //~ bool Ms2(getFlag_("MS2"));
    bool remove_duplicate_features(getFlag_("remove_duplicate_features"));
    
    //-------------------------------------------------------------
    // fetch vocabularies
    //------------------------------------------------------------
    ControlledVocabulary cv;
    cv.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
    cv.loadFromOBO("QC", File::find("/CV/qc-cv.obo"));
 
     QcMLFile qcmlfile;

    //-------------------------------------------------------------
    // MS  aqiusition
    //------------------------------------------------------------
    String base_name = QFileInfo(QString::fromStdString(inputfile_raw)).baseName();

    cout << "Reading mzML file..." << endl;
    MzMLFile mz_data_file;
    MSExperiment<Peak1D> exp;
    MzMLFile().load(inputfile_raw, exp);
    
    //---prep input
    exp.sortSpectra();
    UInt min_mz = std::numeric_limits<UInt>::max();
    UInt max_mz = 0;
    std::map<Size, UInt> mslevelcounts;
    
    qcmlfile.registerRun(base_name,base_name); //TODO use UIDs
    
    //---base MS aquisition qp
    String msaq_ref = base_name + "_msaq";
    QcMLFile::QualityParameter qp;
    qp.id = msaq_ref; ///< Identifier
    qp.cvRef = "QC"; ///< cv reference
    qp.cvAcc = "QC:0000004";
    try
    {
      //~ const ControlledVocabulary::CVTerm& test = cv.getTermByName("MS aquisition result details");
      //~ cout << test.name << test.id << endl;
      const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
      //~ const ControlledVocabulary::CVTerm& term = cv.getTerm("0000004");
      qp.name = term.name; ///< Name
    }
    catch (...)
    {
      qp.name = "mzML file"; ///< Name
    }
    qcmlfile.addRunQualityParameter(base_name, qp);
    
    //---file origin qp
    qp = QcMLFile::QualityParameter();
    qp.name = "mzML file"; ///< Name
    qp.id = base_name + "_run_name"; ///< Identifier
    qp.cvRef = "MS"; ///< cv reference
    qp.cvAcc = "MS:1000577";
    qp.value = base_name;
    qcmlfile.addRunQualityParameter(base_name, qp);
    
    qp = QcMLFile::QualityParameter();
    qp.name = "instrument model"; ///< Name
    qp.id = base_name + "_instrument_name"; ///< Identifier
    qp.cvRef = "MS"; ///< cv reference
    qp.cvAcc = "MS:1000031";
    qp.value = exp.getInstrument().getName();
    qcmlfile.addRunQualityParameter(base_name, qp);    

    qp = QcMLFile::QualityParameter();
    qp.name = "completion time"; ///< Name
    qp.id = base_name + "_date"; ///< Identifier
    qp.cvRef = "MS"; ///< cv reference
    qp.cvAcc = "MS:1000747";
    qp.value = exp.getDateTime().getDate();
    qcmlfile.addRunQualityParameter(base_name, qp);

    //---precursors at
    QcMLFile::Attachment at;
    at.cvRef = "QC"; ///< cv reference
    at.cvAcc = "QC:0000044";
    at.qualityRef = msaq_ref;
    at.id = base_name + "_precursors"; ///< Identifier
    try
    {
      const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
      at.name = term.name; ///< Name
    }
    catch (...)
    {
      at.name = "precursors"; ///< Name
    }

    at.colTypes.push_back("MS:1000894_[sec]"); //RT
    at.colTypes.push_back("MS:1000040"); //MZ
    for (Size i = 0; i < exp.size(); ++i)
    {
      mslevelcounts[exp[i].getMSLevel()]++;
      if (exp[i].getMSLevel() == 2)
      {
        if (exp[i].getPrecursors().front().getMZ() < min_mz)
        {
          min_mz = exp[i].getPrecursors().front().getMZ();
        }
        if (exp[i].getPrecursors().front().getMZ() > max_mz)
        {
          max_mz = exp[i].getPrecursors().front().getMZ();
        }
        std::vector<String> row;
        row.push_back(exp[i].getRT());
        row.push_back(exp[i].getPrecursors().front().getMZ());
        at.tableRows.push_back(row);
      }
    }
    qcmlfile.addRunAttachment(base_name, at);

    //---aquisition results qp
    qp = QcMLFile::QualityParameter();
    qp.cvRef = "QC"; ///< cv reference
    qp.cvAcc = "QC:0000006"; ///< cv accession for "aquisition results"
    qp.id = base_name + "_ms1aquisition"; ///< Identifier
    qp.value = String(mslevelcounts[1]);
    try
    {
      const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
      qp.name = term.name; ///< Name
    }
    catch (...)
    {
      qp.name = "number of ms1 spectra"; ///< Name
    }
    qcmlfile.addRunQualityParameter(base_name, qp);
    

    qp = QcMLFile::QualityParameter();
    qp.cvRef = "QC"; ///< cv reference
    qp.cvAcc = "QC:0000007"; ///< cv accession for "aquisition results"
    qp.id = base_name + "_ms2aquisition"; ///< Identifier
    qp.value = String(mslevelcounts[2]);
    try
    {
      const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
      qp.name = term.name; ///< Name
    }
    catch (...)
    {
      qp.name = "number of ms2 spectra"; ///< Name
    }
    qcmlfile.addRunQualityParameter(base_name, qp);

    qp = QcMLFile::QualityParameter();
    qp.cvRef = "QC"; ///< cv reference
    qp.cvAcc = "QC:0000008"; ///< cv accession for "aquisition results"
    qp.id = base_name + "_Chromaquisition"; ///< Identifier
    qp.value = String(exp.getChromatograms().size());
    try
    {
      const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
      qp.name = term.name; ///< Name
    }
    catch (...)
    {
      qp.name = "number of chromatograms"; ///< Name
    }
    qcmlfile.addRunQualityParameter(base_name, qp);
    
    at = QcMLFile::Attachment();
    at.cvRef = "QC"; ///< cv reference
    at.cvAcc = "QC:0000009";
    at.qualityRef = msaq_ref;
    at.id = base_name + "_mzrange"; ///< Identifier
    try
    {
      const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
      at.name = term.name; ///< Name
    }
    catch (...)
    {
      at.name = "MS MZ aquisition ranges"; ///< Name
    }

    at.colTypes.push_back("QC:0000010"); //MZ
    at.colTypes.push_back("QC:0000011"); //MZ
    std::vector<String> rowmz;
    rowmz.push_back(String(min_mz));
    rowmz.push_back(String(max_mz));
    at.tableRows.push_back(rowmz);
    qcmlfile.addRunAttachment(base_name, at);

    at = QcMLFile::Attachment();
    at.cvRef = "QC"; ///< cv reference
    at.cvAcc = "QC:0000012";
    at.qualityRef = msaq_ref;
    at.id = base_name + "_rtrange"; ///< Identifier
    try
    {
      const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
      at.name = term.name; ///< Name
    }
    catch (...)
    {
      at.name = "MS RT aquisition ranges"; ///< Name
    }

    at.colTypes.push_back("QC:0000013"); //MZ
    at.colTypes.push_back("QC:0000014"); //MZ
    std::vector<String> rowrt;
    rowrt.push_back(String(exp.begin()->getRT()));
    rowrt.push_back(String(exp.getSpectra().back().getRT()));
    at.tableRows.push_back(rowrt);
    qcmlfile.addRunAttachment(base_name, at);
    

    //---ion current stability ( & tic ) qp
    at = QcMLFile::Attachment();
    at.cvRef = "QC"; ///< cv reference
    at.cvAcc = "QC:0000022";
    at.qualityRef = msaq_ref;
    at.id = base_name + "_tics"; ///< Identifier
    try
    {
      const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
      at.name = term.name; ///< Name
    }
    catch (...)
    {
      at.name = "MS TICs"; ///< Name
    }
    
    at.colTypes.push_back("MS:1000894_[sec]");
    at.colTypes.push_back("MS:1000285");
    UInt max = 0;
    Size below_10k = 0;
    for (Size i = 0; i < exp.size(); ++i)
    {
      if (exp[i].getMSLevel() == 1)
      {
        UInt sum = 0;
        for (Size j = 0; j < exp[i].size(); ++j)
        {
          sum += exp[i][j].getIntensity();
        }
        if (sum > max)
        {
          max = sum;
        }
        if (sum < 10000)
        {
          ++below_10k;
        }
        std::vector<String> row;
        row.push_back(exp[i].getRT());
        row.push_back(sum);
        at.tableRows.push_back(row);
      }
    }
    qcmlfile.addRunAttachment(base_name, at);
    

    qp = QcMLFile::QualityParameter();
    qp.id = base_name + "_ticslump"; ///< Identifier
    qp.cvRef = "QC"; ///< cv reference
    qp.cvAcc = "QC:0000023";
    qp.value = String((100 / exp.size()) * below_10k);
    try
    {
      const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
      qp.name = term.name; ///< Name
    }
    catch (...)
    {
      qp.name = "percentage of tic slumps"; ///< Name
    }
    qcmlfile.addRunQualityParameter(base_name, qp);

    
    //-------------------------------------------------------------
    // MS  id
    //------------------------------------------------------------
    if (inputfile_id != "")
    {
      IdXMLFile().load(inputfile_id, prot_ids, pep_ids);
      cerr << "idXML read ended. Found " << pep_ids.size() << " peptide identifications." << endl;

      ProteinIdentification::SearchParameters params = prot_ids[0].getSearchParameters();
      vector<String> var_mods = params.variable_modifications;
      //~ boost::regex re("(?<=[KR])(?=[^P])");
     
      String msid_ref = base_name + "_msid";
      QcMLFile::QualityParameter qp;
      qp.id = msid_ref; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000025";
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
        qp.name = "MS identification result details"; ///< Name
      }
      qcmlfile.addRunQualityParameter(base_name, qp);


      at = QcMLFile::Attachment();
      at.cvRef = "QC"; ///< cv reference
      at.cvAcc = "QC:0000026";
      at.qualityRef = msid_ref;
      at.id = base_name + "_idsetting"; ///< Identifier
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
        at.name = term.name; ///< Name
      }
      catch (...)
      {
        at.name = "MS id settings"; ///< Name
      }
      
      at.colTypes.push_back("MS:1001013"); //MS:1001013 db name  MS:1001016 version  MS:1001020 taxonomy
      at.colTypes.push_back("MS:1001016");
      at.colTypes.push_back("MS:1001020");
      std::vector<String> row;
      row.push_back(String(prot_ids.front().getSearchParameters().db));
      row.push_back(String(prot_ids.front().getSearchParameters().db_version));
      row.push_back(String(prot_ids.front().getSearchParameters().taxonomy));
      at.tableRows.push_back(row);
      qcmlfile.addRunAttachment(base_name, at);


      UInt spectrum_count = 0;
      Size peptide_hit_count = 0;
      UInt runs_count = 0;
      Size protein_hit_count = 0;
      set<String> peptides;
      set<String> proteins;
      Size missedcleavages = 0;
      for (Size i = 0; i < pep_ids.size(); ++i)
      {
        if (!pep_ids[i].empty())
        {
          ++spectrum_count;
          peptide_hit_count += pep_ids[i].getHits().size();
          const vector<PeptideHit>& temp_hits = pep_ids[i].getHits();
          for (Size j = 0; j < temp_hits.size(); ++j)
          {
            peptides.insert(temp_hits[j].getSequence().toString());
          }
        }
      }
      for (set<String>::iterator it = peptides.begin(); it != peptides.end(); ++it)
      {
        for (String::const_iterator st = it->begin(); st != it->end() - 1; ++st)
        {
          if (*st == 'K' || *st == 'R')
          {
            ++missedcleavages;
          }
        }
      }

      for (Size i = 0; i < prot_ids.size(); ++i)
      {
        ++runs_count;
        protein_hit_count += prot_ids[i].getHits().size();
        const vector<ProteinHit>& temp_hits = prot_ids[i].getHits();
        for (Size j = 0; j < temp_hits.size(); ++j)
        {
          proteins.insert(temp_hits[j].getAccession());
        }
      }
      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000037"; ///< cv accession
      qp.id = base_name + "_misscleave"; ///< Identifier
      qp.value = missedcleavages;
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
        qp.name = "total number of missed cleavages"; ///< Name
      }
      qcmlfile.addRunQualityParameter(base_name, qp);

      
      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000032"; ///< cv accession
      qp.id = base_name + "_totprot"; ///< Identifier
      qp.value = protein_hit_count;
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
        qp.name = "total number of identified proteins"; ///< Name
      }
      qcmlfile.addRunQualityParameter(base_name, qp);


      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000033"; ///< cv accession
      qp.id = base_name + "_totuniqprot"; ///< Identifier
      qp.value = String(proteins.size());
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
         qp.name = "total number of uniquely identified proteins"; ///< Name
      }
      qcmlfile.addRunQualityParameter(base_name, qp);


      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000029"; ///< cv accession
      qp.id = base_name + "_psms"; ///< Identifier
      qp.value = String(spectrum_count);
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
         qp.name = "total number of PSM"; ///< Name
      }
      qcmlfile.addRunQualityParameter(base_name, qp);


      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000030"; ///< cv accession
      qp.id = base_name + "_totpeps"; ///< Identifier
      qp.value = String(peptide_hit_count);
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
         qp.name = "total number of identified peptides"; ///< Name
      }
      qcmlfile.addRunQualityParameter(base_name, qp);


      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000031"; ///< cv accession
      qp.id = base_name + "_totuniqpeps"; ///< Identifier
      qp.value = String(peptides.size());
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
         qp.name = "total number of uniquely identified peptides"; ///< Name
      }
      qcmlfile.addRunQualityParameter(base_name, qp);


      at = QcMLFile::Attachment();
      at.cvRef = "QC"; ///< cv reference
      at.cvAcc = "QC:0000038";
      at.qualityRef = msid_ref;
      at.id = base_name + "_massacc"; ///< Identifier
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
        at.name = term.name; ///< Name
      }
      catch (...)
      {
        at.name = "delta ppm tables";
      }
      
      //~ delta ppm QC:0000039 RT MZ uniqueness ProteinID MS:1000885 target/decoy Score PeptideSequence MS:1000889 Annots string Similarity Charge UO:0000219 TheoreticalWeight UO:0000221 Oxidation_(M)
      at.colTypes.push_back("RT");
      at.colTypes.push_back("MZ");
      at.colTypes.push_back("Score");
      at.colTypes.push_back("PeptideSequence");
      at.colTypes.push_back("Charge");
      at.colTypes.push_back("TheoreticalWeight");
      at.colTypes.push_back("delta_ppm");
      for (UInt w = 0; w < var_mods.size(); ++w)
      {
        at.colTypes.push_back(String(var_mods[w]).substitute(' ', '_'));
      }

      std::vector<double> deltas;
      //~ prot_ids[0].getSearchParameters();
      for (vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
      {
        if (it->getHits().size() > 0)
        {
          std::vector<String> row;
          row.push_back(it->getRT());
          row.push_back(it->getMZ());
          PeptideHit tmp = it->getHits().front(); //TODO depends on score & sort
          vector<UInt> pep_mods;
          for (UInt w = 0; w < var_mods.size(); ++w)
          {
            pep_mods.push_back(0);
          }
          for (AASequence::ConstIterator z =  tmp.getSequence().begin(); z != tmp.getSequence().end(); ++z)
          {
            Residue res = *z;
            String temp;
            if (res.getModification().size() > 0 && res.getModification() != "Carbamidomethyl")
            {
              temp = res.getModification() + " (" + res.getOneLetterCode()  + ")";
              //cout<<res.getModification()<<endl;
              for (UInt w = 0; w < var_mods.size(); ++w)
              {
                if (temp == var_mods[w])
                {
                  //cout<<temp;
                  pep_mods[w] += 1;
                }
              }
            }
          }
          row.push_back(tmp.getScore());
          row.push_back(tmp.getSequence().toString().removeWhitespaces());
          row.push_back(tmp.getCharge());
          row.push_back(String((tmp.getSequence().getMonoWeight() + tmp.getCharge() * Constants::PROTON_MASS_U) / tmp.getCharge()));
          double dppm = /* std::abs */ (getMassDifference(((tmp.getSequence().getMonoWeight() + tmp.getCharge() * Constants::PROTON_MASS_U) / tmp.getCharge()), it->getMZ(), true));
          row.push_back(String(dppm));
          deltas.push_back(dppm);
          for (UInt w = 0; w < var_mods.size(); ++w)
          {
            row.push_back(pep_mods[w]);
          }
          at.tableRows.push_back(row);
        }
      }
      qcmlfile.addRunAttachment(base_name, at);
      

      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000040"; ///< cv accession
      qp.id = base_name + "_mean_delta"; ///< Identifier
      qp.value = String(OpenMS::Math::mean(deltas.begin(), deltas.end()));
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
         qp.name = "mean delta ppm"; ///< Name
      }
      qcmlfile.addRunQualityParameter(base_name, qp);


      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000041"; ///< cv accession
      qp.id = base_name + "_median_delta"; ///< Identifier
      qp.value = String(OpenMS::Math::median(deltas.begin(), deltas.end(), false));
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
         qp.name = "median delta ppm"; ///< Name
      }
      qcmlfile.addRunQualityParameter(base_name, qp);


      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000035"; ///< cv accession
      qp.id = base_name + "_ratio_id"; ///< Identifier
      qp.value = String(double(pep_ids.size()) / double(mslevelcounts[2]));
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
         qp.name = "id ratio"; ///< Name
      }
      qcmlfile.addRunQualityParameter(base_name, qp);
    }

    //-------------------------------------------------------------
    // MS quantitation
    //------------------------------------------------------------
    FeatureMap map;
    String msqu_ref = base_name + "_msqu";
    if (inputfile_feature != "")
    {
      FeatureXMLFile f;
      f.load(inputfile_feature, map);

      cout << "Read featureXML file..." << endl;

      //~ UInt fiter = 0;
      map.sortByRT();
      map.updateRanges();

      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000045"; ///< cv accession
      qp.id = msqu_ref; ///< Identifier
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
         qp.name = "MS quantification result details"; ///< Name
      }
      qcmlfile.addRunQualityParameter(base_name, qp);
      
      qp = QcMLFile::QualityParameter();
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:0000046"; ///< cv accession
      qp.id = base_name + "_feature_count"; ///< Identifier
      qp.value = String(map.size());
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(qp.cvAcc);
        qp.name = term.name; ///< Name
      }
      catch (...)
      {
         qp.name = "number of features"; ///< Name
      }
      qcmlfile.addRunQualityParameter(base_name, qp);      
    }

    if (inputfile_feature != "" && !remove_duplicate_features)
    {
      
      QcMLFile::Attachment at;
      at = QcMLFile::Attachment();
      at.cvRef = "QC"; ///< cv reference
      at.cvAcc = "QC:0000047";
      at.qualityRef = msqu_ref;
      at.id = base_name + "_features"; ///< Identifier
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
        at.name = term.name; ///< Name
      }
      catch (...)
      {
        at.name = "features"; ///< Name
      }
      
      at.colTypes.push_back("MZ");
      at.colTypes.push_back("RT");
      at.colTypes.push_back("Intensity");
      at.colTypes.push_back("Charge");
      at.colTypes.push_back("Quality");
      at.colTypes.push_back("FWHM");
      at.colTypes.push_back("IDs");
      UInt fiter = 0;
      map.sortByRT();
      //ofstream out(outputfile_name.c_str());
      while (fiter < map.size())
      {
        std::vector<String> row;
        row.push_back(map[fiter].getMZ());
        row.push_back(map[fiter].getRT());
        row.push_back(map[fiter].getIntensity());
        row.push_back(map[fiter].getCharge());
        row.push_back(map[fiter].getOverallQuality());
        row.push_back(map[fiter].getWidth());
        row.push_back(map[fiter].getPeptideIdentifications().size());
        fiter++;
        at.tableRows.push_back(row);
      }     
      qcmlfile.addRunAttachment(base_name, at);
    }
    else if (inputfile_feature != "" && remove_duplicate_features)
    {
      QcMLFile::Attachment at;
      at = QcMLFile::Attachment();
      at.cvRef = "QC"; ///< cv reference
      at.cvAcc = "QC:0000047";
      at.qualityRef = msqu_ref;
      at.id = base_name + "_features"; ///< Identifier
      try
      {
        const ControlledVocabulary::CVTerm& term = cv.getTerm(at.cvAcc);
        at.name = term.name; ///< Name
      }
      catch (...)
      {
        at.name = "features"; ///< Name
      }
      
      at.colTypes.push_back("MZ");
      at.colTypes.push_back("RT");
      at.colTypes.push_back("Intensity");
      at.colTypes.push_back("Charge");
      FeatureMap map, map_out;
      FeatureXMLFile f;
      f.load(inputfile_feature, map);
      UInt fiter = 0;
      map.sortByRT();
      while (fiter < map.size())
      {
        FeatureMap map_tmp;
        for (UInt k = fiter; k <= map.size(); ++k)
        {
          if (abs(map[fiter].getRT() - map[k].getRT()) < 0.1)
          {
            //~ cout << fiter << endl;
            map_tmp.push_back(map[k]);
          }
          else
          {
            fiter = k;
            break;
          }
        }
        map_tmp.sortByMZ();
        UInt retif = 1;
        map_out.push_back(map_tmp[0]);
        while (retif < map_tmp.size())
        {
          if (abs(map_tmp[retif].getMZ() - map_tmp[retif - 1].getMZ()) > 0.01)
          {
            cout << "equal RT, but mass different" << endl;
            map_out.push_back(map_tmp[retif]);
          }
          retif++;
        }
      }
      qcmlfile.addRunAttachment(base_name, at);
    }
    if (inputfile_consensus != "")
    {
      cout << "Reading consensusXML file..." << endl;
      ConsensusXMLFile f;
      ConsensusMap map;
      f.load(inputfile_consensus, map);
      //~ String CONSENSUS_NAME = "_consensus.tsv";
      //~ String combined_out = outputfile_name + CONSENSUS_NAME;
      //~ ofstream out(combined_out.c_str());

      at = QcMLFile::Attachment();
      qp.name = "consensuspoints"; ///< Name
      //~ qp.id = base_name + "_consensuses"; ///< Identifier
      qp.cvRef = "QC"; ///< cv reference
      qp.cvAcc = "QC:xxxxxxxx"; ///< cv accession "featuremapper results"

      at.colTypes.push_back("Native_spectrum_ID");
      at.colTypes.push_back("DECON_RT_(sec)");
      at.colTypes.push_back("DECON_MZ_(Th)");
      at.colTypes.push_back("DECON_Intensity");
      at.colTypes.push_back("Feature_RT_(sec)");
      at.colTypes.push_back("Feature_MZ_(Th)");
      at.colTypes.push_back("Feature_Intensity");
      at.colTypes.push_back("Feature_Charge");
      for (ConsensusMap::const_iterator cmit = map.begin(); cmit != map.end(); ++cmit)
      {
        const ConsensusFeature& CF = *cmit;
        for (ConsensusFeature::const_iterator cfit = CF.begin(); cfit != CF.end(); ++cfit)
        {
          std::vector<String> row;
          FeatureHandle FH = *cfit;
          row.push_back(CF.getMetaValue("spectrum_native_id"));
          row.push_back(CF.getRT()); row.push_back(CF.getMZ());
          row.push_back(CF.getIntensity());
          row.push_back(FH.getRT());
          row.push_back(FH.getMZ());
          row.push_back(FH.getCharge());
          at.tableRows.push_back(row);
        }
      }
      qcmlfile.addRunAttachment(base_name, at);
    }
    
    
    //-------------------------------------------------------------
    // finalize
    //------------------------------------------------------------
    qcmlfile.store(outputfile_name);
    return EXECUTION_OK;
  }

};

#pragma clang diagnostic pop

int main(int argc, const char** argv)
{
  TOPPQCCalculator tool;
  return tool.main(argc, argv);
}

/// @endcond
