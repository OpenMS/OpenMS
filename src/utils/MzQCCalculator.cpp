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
// $Maintainer: Axel Walter $
// $Author: Mathias Walzer, Sven Nahnsen $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzQCFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/SYSTEM/File.h>



#include <QFileInfo>
//~ #include <boost/regex.hpp>

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MzQCCalculator MzQCCalculator

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
    
    - @p id produces quality parameter values for the identification file; this file should contain either only the final psm to each spectrum (1 PeptideHit per identified spectrum) or have the PeptideHits sorted to 'best' first, where 'best' depends on the use case.
    - @p feature produces quality parameter values for the feature file; this file can be either mapped or unmapped, the latter reulting in less metrics available.
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

class MzQCCalculator :
  public TOPPBase
{
public:

  MzQCCalculator() :
    TOPPBase("QCCalculator", 
      "Calculates basic quality parameters from MS experiments and subsequent analysis data as identification or feature detection.", 
      false, 
      {{ "Walzer M, Pernas LE, Nasso S, Bittremieux W, Nahnsen S, Kelchtermans P,  Martens, L", 
         "qcML: An Exchange Format for Quality Control Metrics from Mass Spectrometry Experiments", 
         "Molecular & Cellular Proteomics 2014; 13(8)" , "10.1074/mcp.M113.035907"
      }})
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "raw data input file (this is relevant if you want to look at MS1, MS2 and precursor peak information)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Your mzQC file.");
    setValidFormats_("out", ListUtils::create<String>("qcML"));
    registerStringOption_("label", "<label>", "", "unique name for the run that can be used in a figure label", true);
    registerStringOption_("name", "<contactName>", "", "name of the person creating this mzQC file", false);
    registerStringOption_("address", "<contactAddress>", "", "contact address (mail/e-mail or phone)", false);
    registerStringOption_("description", "<description>", "", "description and comments about the mzQC file contents", false);
    registerInputFile_("id", "<file>", "", "Input idXML file containing the identifications. Your identifications will be exported in an easy-to-read format", false);
    setValidFormats_("id", ListUtils::create<String>("idXML"));
    registerInputFile_("feature", "<file>", "", "feature input file (this is relevant for most QC issues)", false);
    setValidFormats_("feature", ListUtils::create<String>("featureXML"));
    registerInputFile_("consensus", "<file>", "", "consensus input file (this is only used for charge state deconvoluted output. Use the consensusXML output form the DeCharger)", false);
    setValidFormats_("consensus", ListUtils::create<String>("consensusXML"));
    registerFlag_("remove_duplicate_features", "This flag should be set, if you work with a set of merged features.");
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

  float calculateSNmedian (MSSpectrum& spec, bool norm = true)
  {
    if (spec.size() == 0) return 0;
    float median = 0;
    float maxi = 0;
    spec.sortByIntensity();
    
    if (spec.size() % 2 == 0)
    {
      median = (spec[spec.size() / 2 - 1].getIntensity() + spec[spec.size() / 2].getIntensity()) / 2;
    }
    else
    {
      median = spec[spec.size() / 2].getIntensity();
    }
    maxi = spec.back().getIntensity();
    if (!norm)
    {
      float sn_by_max2median = maxi / median;
      return sn_by_max2median;
    }

    float sign_int= 0;
    float nois_int = 0;
    size_t sign_cnt= 0;
    size_t nois_cnt = 0;
    for (MSSpectrum::const_iterator pt = spec.begin(); pt != spec.end(); ++pt)
    {
      if (pt->getIntensity() <= median)
      {
        ++nois_cnt;
        nois_int += pt->getIntensity();
      }
      else
      {
        ++sign_cnt;
        sign_int += pt->getIntensity();
      }
    }
    float sn_by_max2median_norm = (sign_int / sign_cnt) / (nois_int / nois_cnt);

    return sn_by_max2median_norm;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_feature = getStringOption_("feature");
    String inputfile_consensus = getStringOption_("consensus");
    String inputFileName = getStringOption_("in");
    String outputFileName = getStringOption_("out");
    String contactName = getStringOption_("name");
    String contactAddress = getStringOption_("address");
    String description = getStringOption_("description");
    String label = getStringOption_("label");
    
    // bool remove_duplicate_features(getFlag_("remove_duplicate_features"));
    
    //-------------------------------------------------------------
    // fetch vocabularies
    //------------------------------------------------------------
    ControlledVocabulary cv;
    cv.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
    cv.loadFromOBO("QC", File::find("/CV/qc-cv.obo"));
 
    MzQCFile qcfile;

    cout << "Reading mzML file..." << endl;
    PeakMap exp;
    MzMLFile().load(inputFileName, exp);
    exp.sortSpectra();
    // UInt min_mz = std::numeric_limits<UInt>::max();
    //  UInt max_mz = 0;
    std::map<Size, UInt> mslevelcounts;

    qcfile.store(inputFileName, outputFileName, exp, cv, contactName, contactAddress, description, label);
    return EXECUTION_OK;
  }

};

#pragma clang diagnostic pop

int main(int argc, const char** argv)
{
  MzQCCalculator tool;
  return tool.main(argc, argv);
}

/// @endcond
