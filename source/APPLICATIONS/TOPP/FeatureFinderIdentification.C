// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
        @page TOPP_FeatureFinderIdentification FeatureFinderIdentification

        @brief Detects features in MS1 data based on peptide identifications.

        <CENTER>
        <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderIdentification \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_MapAlignerIdentification</td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
        </tr>
        </table>
        </CENTER>

        This tool uses algorithms for targeted data analysis from the OpenSWATH pipeline.

        <B>The command line parameters of this tool are:</B>
        @verbinclude TOPP_FeatureFinderIdentification.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderIdentification :
  public TOPPBase
{
public:
  TOPPFeatureFinderIdentification() :
    TOPPBase("FeatureFinderIdentification", "Detects features in MS1 data based on peptide identifications.", false)
  {
  }

protected:

  typedef MSExperiment<Peak1D> PeakMap;

  // mapping: charge -> iterator to peptide
  typedef Map<Int, vector<vector<PeptideIdentification>::iterator> > ChargeMap;
  // mapping: sequence -> charge -> iterator to peptide
  typedef Map<AASequence, ChargeMap> PeptideMap;

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file (LC-MS raw data)");
    setValidFormats_("in", StringList::create("mzML"));
    registerInputFile_("id", "<file>", "", 
                       "input file (peptide identifications)");
    setValidFormats_("id", StringList::create("idXML"));
    registerOutputFile_("out", "<file>", "", "output file (features)");
    setValidFormats_("out", StringList::create("featureXML"));
    registerOutputFile_("lib_out","<file>", "", "output file (library)", false);
    setValidFormats_("lib_out", StringList::create("traML"));
    registerOutputFile_("chrom_out","<file>", "", "output file (chromatograms)",
                        false);
    setValidFormats_("chrom_out", StringList::create("mzML"));
    registerOutputFile_("trafo_out","<file>", "", 
                        "output file (RT transformation)", false);
    setValidFormats_("trafo_out", StringList::create("trafoXML"));

    addEmptyLine_();
    registerDoubleOption_("isotope_pmin", "<value>", 0.01, "Minimum probability for an isotope to be included in the assay for a peptide.", false);
    setMinFloat_("isotope_pmin", 0);
    setMaxFloat_("isotope_pmin", 1);
    registerDoubleOption_("rt_window", "<value>", 180, "RT window size (in sec.) for chromatogram extraction.", false);
    setMinFloat_("rt_window", 0);
    registerDoubleOption_("mz_window", "<value>", 0.05, "m/z window size (in Th) for chromatogram extraction.", false);
    setMinFloat_("mz_window", 0);

    // addEmptyLine_();
    // registerSubsection_("algorithm", "Algorithm parameters section");
  }


  // Param getSubsectionDefaults_(const String& /*section*/) const
  // {
  //   Param combined;
  //   return combined;
  // }


  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String id = getStringOption_("id");
    String out = getStringOption_("out");
    String lib_out = getStringOption_("lib_out");
    String chrom_out = getStringOption_("chrom_out");
    String trafo_out = getStringOption_("trafo_out");
    DoubleReal isotope_pmin = getDoubleOption_("isotope_pmin");
    DoubleReal rt_window = getDoubleOption_("rt_window");
    DoubleReal mz_window = getDoubleOption_("mz_window");

    //-------------------------------------------------------------
    // load input
    //-------------------------------------------------------------
    LOG_INFO << "Loading input data..." << endl;
    MzMLFile mzml;
    mzml.setLogType(log_type_);
    mzml.getOptions().addMSLevel(1);
    PeakMap ms_data;
    mzml.load(in, ms_data);

    // RT transformation to range 0-1:
    ms_data.updateRanges();
    DoubleReal min_rt = ms_data.getMinRT(), max_rt = ms_data.getMaxRT();
    TransformationDescription trafo;
    TransformationDescription::DataPoints points;
    points.push_back(make_pair(min_rt, 0.0));
    points.push_back(make_pair(max_rt, 1.0));
    trafo.setDataPoints(points);
    trafo.fitModel("linear");
    if (!trafo_out.empty())
    {
      TransformationXMLFile().store(trafo_out, trafo);
    }
    
    vector<PeptideIdentification> peptides;
    vector<ProteinIdentification> proteins;
    IdXMLFile().load(id, proteins, peptides);

    //-------------------------------------------------------------
    // prepare peptide map
    //-------------------------------------------------------------
    LOG_INFO << "Preparing mapping of peptide data..." << endl;
    PeptideMap peptide_map;
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin(); 
         pep_it != peptides.end(); ++pep_it)
    {
      if (pep_it->getHits().empty()) continue;
      pep_it->sort();
      PeptideHit& hit = pep_it->getHits()[0];
      peptide_map[hit.getSequence()][hit.getCharge()].push_back(pep_it);
    }

    //-------------------------------------------------------------
    // create assay library from peptides
    //-------------------------------------------------------------
    LOG_INFO << "Creating assay library..." << endl;
    TargetedExperiment library;
    set<String> protein_accessions;

    for (PeptideMap::iterator pm_it = peptide_map.begin(); 
         pm_it != peptide_map.end(); ++pm_it)
    {
      const AASequence& seq = pm_it->first;
      // LOG_DEBUG << "Peptide: " << seq.toString() << endl;

      // keep track of protein accessions:
      const PeptideHit& hit = pm_it->second.begin()->second[0]->getHits()[0];
      vector<String> current_accessions = hit.getProteinAccessions();
      // missing protein accession would crash OpenSwath algorithms:
      if (current_accessions.empty())
      {
        current_accessions.push_back("not_available");
      }
      protein_accessions.insert(current_accessions.begin(), 
                                current_accessions.end());

      // get isotope distribution for peptide:
      IsotopeDistribution iso_dist = 
        seq.getFormula(Residue::Full, 0).getIsotopeDistribution(10);
      iso_dist.trimLeft(isotope_pmin);
      iso_dist.trimRight(isotope_pmin);
      iso_dist.renormalize();

      // go through different charge states:
      for (ChargeMap::iterator cm_it = pm_it->second.begin(); 
           cm_it != pm_it->second.end(); ++cm_it)
      {
        Int charge = cm_it->first;
        DoubleReal mz = seq.getMonoWeight(Residue::Full, charge) / charge;

        // get median RT and normalize it:
        DoubleList rts;
        for (vector<vector<PeptideIdentification>::iterator>::iterator pi_it = 
               cm_it->second.begin(); pi_it != cm_it->second.end(); ++pi_it)
        {
          rts << (*pi_it)->getMetaValue("RT");
        }
        DoubleReal median_rt = Math::median(rts.begin(), rts.end());
        CVTerm rt_term;
        rt_term.setCVIdentifierRef("MS");
        rt_term.setAccession("MS:1000896");
        rt_term.setName("normalized retention time");
        rt_term.setValue(trafo.apply(median_rt));

        // create assay for current peptide and charge state:
        TargetedExperiment::Peptide peptide;
        peptide.sequence = seq.toString();
        peptide.id = peptide.sequence + "/" + String(charge);
        peptide.protein_refs = current_accessions;
        peptide.setChargeState(charge);
        TargetedExperiment::RetentionTime rt;
        rt.addCVTerm(rt_term);
        peptide.rts.push_back(rt);
        library.addPeptide(peptide);
        
        // go through different isotopes, add transitions:
        Size counter = 0;
        for (IsotopeDistribution::iterator iso_it = iso_dist.begin();
             iso_it != iso_dist.end(); ++iso_it, ++counter)
        {
          String annotation = "i" + String(counter);
          String transition_name = peptide.id + "_" + annotation;
          
          ReactionMonitoringTransition transition;
          transition.setNativeID(transition_name);
          transition.setPrecursorMZ(mz);
          transition.setProductMZ(mz + float(counter) / charge);
          transition.setLibraryIntensity(iso_it->second * 100);
          transition.setMetaValue("annotation", annotation);
          transition.setPeptideRef(peptide.id);
          library.addTransition(transition);
        }
      }      
    }
    // add protein references:
    for (set<String>::iterator acc_it = protein_accessions.begin();
         acc_it != protein_accessions.end(); ++acc_it)
    {
      TargetedExperiment::Protein protein;
      protein.id = *acc_it;
      library.addProtein(protein);
    }

    if (!lib_out.empty())
    {
      TraMLFile().store(lib_out, library);
    }

    //-------------------------------------------------------------
    // extract chromatograms
    //-------------------------------------------------------------
    LOG_INFO << "Extracting chromatograms..." << endl;
    ChromatogramExtractor extractor;
    PeakMap chrom_data;
    extractor.setLogType(log_type_);
    extractor.extractChromatograms(ms_data, chrom_data, library, mz_window,
                                   false, trafo, rt_window / 2.0, "tophat");
    if (!chrom_out.empty())
    {
      MzMLFile().store(chrom_out, chrom_data);
    }

    //-------------------------------------------------------------
    // find chromatographic peaks
    //-------------------------------------------------------------
    LOG_INFO << "Finding chromatographic peaks..." << endl;
    FeatureMap<> features;
    PeakMap dummy;
    MRMFeatureFinderScoring mrm_finder;
    Param params = mrm_finder.getParameters();
    params.setValue("stop_report_after_feature", 1);
    params.setValue("TransitionGroupPicker:PeakPickerMRM:use_gauss", "false");
    params.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
    params.setValue("TransitionGroupPicker:PeakPickerMRM:method", "corrected");
    mrm_finder.setParameters(params);
    mrm_finder.setLogType(log_type_);
    mrm_finder.setStrictFlag(false);
    mrm_finder.pickExperiment(chrom_data, features, library, trafo, dummy);

    //-------------------------------------------------------------
    // fill in missing feature data
    //-------------------------------------------------------------
    LOG_INFO << "Adapting feature data..." << endl;
    for (FeatureMap<>::Iterator feat_it = features.begin(); 
         feat_it != features.end(); ++feat_it)
    {
      feat_it->setMZ(feat_it->getMetaValue("PrecursorMZ"));
      feat_it->setCharge(feat_it->getPeptideIdentifications()[0].getHits()[0].
                         getCharge());
      DoubleReal rt_min = feat_it->getMetaValue("leftWidth");
      DoubleReal rt_max = feat_it->getMetaValue("rightWidth");
      if (feat_it->getConvexHulls().empty()) // add hulls for mass traces
      {
        for (vector<Feature>::iterator sub_it = 
               feat_it->getSubordinates().begin(); sub_it !=
               feat_it->getSubordinates().end(); ++sub_it)
        {
          ConvexHull2D hull;
          hull.addPoint(DPosition<2>(rt_min, sub_it->getMZ() - mz_window / 2));
          hull.addPoint(DPosition<2>(rt_min, sub_it->getMZ() + mz_window / 2));
          hull.addPoint(DPosition<2>(rt_max, sub_it->getMZ() - mz_window / 2));
          hull.addPoint(DPosition<2>(rt_max, sub_it->getMZ() + mz_window / 2));
          feat_it->getConvexHulls().push_back(hull);
        }
      }
    }

    //-------------------------------------------------------------
    // write output
    //-------------------------------------------------------------
    LOG_INFO << "Writing results..." << endl;
    features.ensureUniqueId();
    addDataProcessing_(features, 
                       getProcessingInfo_(DataProcessing::QUANTITATION));
    FeatureXMLFile().store(out, features);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderIdentification tool;
  return tool.main(argc, argv);
}

/// @endcond
