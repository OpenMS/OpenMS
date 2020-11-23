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
// $Maintainer: Chris Bielow $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>
#include <OpenMS/QC/Contaminants.h>
#include <OpenMS/QC/FragmentMassError.h>
#include <OpenMS/QC/FWHM.h>
#include <OpenMS/QC/MissedCleavages.h>
#include <OpenMS/QC/Ms2IdentificationRate.h>
#include <OpenMS/QC/MzCalibration.h>
#include <OpenMS/QC/PeptideMass.h>
#include <OpenMS/QC/RTAlignment.h>
#include <OpenMS/QC/TIC.h>
#include <OpenMS/QC/Ms2SpectrumStats.h>
#include <cstdio>

#include <map>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_QualityControl QualityControl

@brief Generates an mzTab file from various sources of a pipeline (mainly a ConsensusXML) which can be used for QC plots (e.g. via the R package 'PTXQC').

<CENTER>
<table>
<tr>
<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ QualityControl \f$ \longrightarrow \f$</td>
<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeledKD (or FLs; for consensusXML)</td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=3> <a href="https://github.com/cbielow/PTXQC/" target="_blank">PTX-QC</a> </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper (for featureXMLs)</td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_InternalCalibration </td>
</tr>
</table>
</CENTER>

See @ref TOPP_example_qualitycontrol for details.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_QualityControl.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_QualityControl.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPQualityControl : public TOPPBase
{
public:
  TOPPQualityControl()
   : TOPPBase("QualityControl", "Computes various QC metrics from many possible input files (only the consensusXML is required). The more optional files you provide, the more metrics you get.", true)
  {
  }
protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_cm", "<file>", "", "ConsensusXML input, generated by FeatureLinker.", true);
    setValidFormats_("in_cm", {"consensusXML"});
    registerInputFileList_("in_raw", "<files>", {}, "MzML input (after InternalCalibration, if available)", false);
    setValidFormats_("in_raw", {"mzML"});
    registerInputFileList_("in_postFDR", "<files>", {}, "FeatureXMLs after FDR filtering", false);
    setValidFormats_("in_postFDR", {"featureXML"});
    registerOutputFile_("out", "<file>", "", "Output mzTab with QC information", false);
    setValidFormats_("out", { "mzTab" });
    registerOutputFile_("out_cm", "<file>", "", "ConsensusXML with QC information (as metavalues)", false);
    setValidFormats_("out_cm", { "consensusXML" });
    registerOutputFileList_("out_feat", "<files>", {}, "FeatureXMLs with QC information (as metavalues)", false);
    setValidFormats_("out_feat", { "featureXML" });
    registerTOPPSubsection_("FragmentMassError", "Fragment Mass Error settings");
    registerStringOption_("FragmentMassError:unit", "<unit>", "auto", "Unit for mass tolerance. 'auto' uses information from FeatureXML", false);
    setValidStrings_("FragmentMassError:unit", std::vector<String>(FragmentMassError::names_of_toleranceUnit, FragmentMassError::names_of_toleranceUnit + (int)FragmentMassError::ToleranceUnit::SIZE_OF_TOLERANCEUNIT));
    registerDoubleOption_("FragmentMassError:tolerance", "<double>", 20, "m/z search window for matching peaks in two spectra", false);
    registerInputFile_("in_contaminants", "<file>", "", "Proteins considered contaminants", false);
    setValidFormats_("in_contaminants", {"fasta"});
    registerInputFileList_("in_trafo", "<file>", {}, "trafoXMLs from MapAligners", false);
    setValidFormats_("in_trafo", {"trafoXML"});
    registerTOPPSubsection_("MS2_id_rate", "MS2 ID Rate settings");
    registerFlag_("MS2_id_rate:assume_all_target", "Forces the metric to run even if target/decoy annotation is missing (accepts all pep_ids as target hits).", false);
    //TODO get ProteinQuantifier output for PRT section
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    //
    // Read input, check for same length and get that length
    QCBase::Status status;
    UInt64 number_exps(0);
    StringList in_raw = updateFileStatus_(status, number_exps, "in_raw", QCBase::Requires::RAWMZML);
    StringList in_postFDR = updateFileStatus_(status, number_exps, "in_postFDR", QCBase::Requires::POSTFDRFEAT);
    StringList in_trafo = updateFileStatus_(status, number_exps, "in_trafo", QCBase::Requires::TRAFOALIGN);

    // load databases and other single file inputs
    String in_contaminants = getStringOption_("in_contaminants");
    vector<FASTAFile::FASTAEntry> contaminants;
    if (!in_contaminants.empty())
    {
      FASTAFile::load(in_contaminants, contaminants);
      status |= QCBase::Requires::CONTAMINANTS;
    }
    ConsensusMap cmap;
    String in_cm = getStringOption_("in_cm");
    ConsensusXMLFile().load(in_cm, cmap);

    std::vector<FeatureMap> fmaps;
    if (in_postFDR.empty())
    {
      status |= QCBase::Requires::POSTFDRFEAT;
      fmaps = cmap.split(ConsensusMap::SplitMeta::COPY_ALL);
      bool is_labeled_cmap = QCBase::isLabeledExperiment(cmap);
      if (is_labeled_cmap) // for labeled input (e.g. iTRAQ/TMT/SILAC)
      {
        OPENMS_LOG_INFO << "Labeled data detected!" << std::endl;
        if (number_exps != 1) // no features given, but >1 trafos...
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("More than one mzML or TrafoXML were given, but this is not supported in 'labeled' mode."));
        }
        // number_exps can remain 1, since we only need to annotate the first FMap with metavalues (the others only have exact copies)
        // ...
      }
      else // unlabeled == LFQ mode
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        // currently missing:
        // - invert RT of all features+their PepIDs to allow RTmetric to work (if TrafoXMLs are provided) -- or even better: delegate this to the RTMetric
        // - the SearchParameters are currently taken from the first ProteinIdentificaion of the FMaps...
        //    however, during splitting, all ProtID's from the CMap are blindly copied to all FMaps (it should only pick the correct one)...

        OPENMS_LOG_INFO << "Unlabeled data detected in ConsensusXML detected! Data will be extracted from there. If you can, provide the FeatureXML files for potentially more metrics." << std::endl;
        if (number_exps != fmaps.size())
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        String("Number of Maps in the ConsensusMap (") + fmaps.size() +
                                        ") does not match length of -in_raw or -in_trafo (" + number_exps + ").");
        }
      }
    }
    FeatureMap* fmap;

    //-------------------------------------------------------------
    // prot/pepID-identifier -->  ms-run-path
    //-------------------------------------------------------------
    ProteinIdentification::Mapping mp_c(cmap.getProteinIdentifications());

    //-------------------------------------------------------------
    // Build a PepID Map to later find the corresponding PepID in the CMap
    //-------------------------------------------------------------
    multimap<String, PeptideIdentification*> customID_to_cpepID; // multimap is required because a PepID could be duplicated by IDMapper and appear >=1 in a featureMap
    for (Size i = 0; i < cmap.size(); ++i)
    {
      fillConsensusPepIDMap_(cmap[i].getPeptideIdentifications(), mp_c.identifier_to_msrunpath, customID_to_cpepID);
      // connect CF (stored in PEP section) with its peptides (stored in PSM section) ... they might get separated later by IDConflictResolverAlgorithm
      cmap[i].setMetaValue("cf_id", i);
      for (auto& pep_id : cmap[i].getPeptideIdentifications()) pep_id.setMetaValue("cf_id", i);
    }
    fillConsensusPepIDMap_(cmap.getUnassignedPeptideIdentifications(), mp_c.identifier_to_msrunpath, customID_to_cpepID);
    for (auto& pep_id : cmap.getUnassignedPeptideIdentifications()) pep_id.setMetaValue("cf_id", -1);


    // check flags
    bool all_target_flag = getFlag_("MS2_id_rate:assume_all_target");
    double tolerance_value = getDoubleOption_("FragmentMassError:tolerance");

    auto it = std::find(FragmentMassError::names_of_toleranceUnit, FragmentMassError::names_of_toleranceUnit + (int)FragmentMassError::ToleranceUnit::SIZE_OF_TOLERANCEUNIT, getStringOption_("FragmentMassError:unit"));
    auto idx = std::distance(FragmentMassError::names_of_toleranceUnit, it);
    auto tolerance_unit = FragmentMassError::ToleranceUnit(idx);


    // Instantiate the QC metrics
    Contaminants qc_contaminants;
    FragmentMassError qc_frag_mass_err;
    FWHM qc_fwhm;
    MissedCleavages qc_missed_cleavages;
    Ms2IdentificationRate qc_ms2ir;
    MzCalibration qc_mz_calibration;
    RTAlignment qc_rt_alignment;
    PeptideMass qc_pepmass;
    TIC qc_tic;
    Ms2SpectrumStats qc_ms2stats;
    MzMLFile mzml_file;
    PeakMap exp;
    QCBase::SpectraMap spec_map;

    // Loop through featuremaps...
    vector<PeptideIdentification> all_new_upep_ids;
    for (Size i = 0; i < number_exps; ++i)
    {
      //-------------------------------------------------------------
      // reading input
      //-------------------------------------------------------------
      if (i < in_raw.size())
      { // we either have 'n' or 1 mzML ... use the correct one in each iteration
        mzml_file.load(in_raw[i], exp);
        spec_map.calculateMap(exp);
      }

      ProteinIdentification::Mapping mp_f;
      FeatureXMLFile fxml_file;
      FeatureMap fmap_local;
      if (!in_postFDR.empty())
      {
        fxml_file.load(in_postFDR[i], fmap_local);
        fmap = &fmap_local;
      }
      else
      {
        fmap = &(fmaps[i]);
      }
      mp_f.create(fmap->getProteinIdentifications());

      TransformationXMLFile trafo_file;
      TransformationDescription trafo_descr;
      if (!in_trafo.empty())
      {
        trafo_file.load(in_trafo[i], trafo_descr);
      }
      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------

      if (qc_contaminants.isRunnable(status))
      {
        qc_contaminants.compute(*fmap, contaminants);
      }

      if (qc_frag_mass_err.isRunnable(status))
      {
        qc_frag_mass_err.compute(*fmap, exp, spec_map, tolerance_unit, tolerance_value);
      }

      if (qc_ms2ir.isRunnable(status))
      {
        qc_ms2ir.compute(*fmap, exp, all_target_flag);
      }

      if (qc_mz_calibration.isRunnable(status))
      {
        qc_mz_calibration.compute(*fmap, exp, spec_map);
      }
      
      // after qc_mz_calibration, because it calculates 'mass' metavalue
      if (qc_missed_cleavages.isRunnable(status))
      {
        qc_missed_cleavages.compute(*fmap);
      }

      if (qc_rt_alignment.isRunnable(status))
      { // add metavalues rt_raw & rt_align to all PepIDs
        qc_rt_alignment.compute(*fmap, trafo_descr);
      }

      if (qc_fwhm.isRunnable(status))
      {
        qc_fwhm.compute(*fmap);
      }

      if (qc_pepmass.isRunnable(status))
      {
        qc_pepmass.compute(*fmap);
      }

      if (qc_tic.isRunnable(status))
      {
        qc_tic.compute(exp);
      }

      if (qc_ms2stats.isRunnable(status))
      {
        // copies FWHM metavalue to PepIDs as well
        vector<PeptideIdentification> new_upep_ids = qc_ms2stats.compute(exp, *fmap, spec_map);
        // use identifier of CMap for just calculated pepIDs (via common MS-run-path)
        const auto& f_runpath = mp_f.runpath_to_identifier.begin()->first; // just get any runpath from fmap
        const auto ptr_cmap = mp_c.runpath_to_identifier.find(f_runpath);
        if (ptr_cmap == mp_c.runpath_to_identifier.end())
        {
          OPENMS_LOG_ERROR << "FeatureXML (MS run '" << ListUtils::concatenate(f_runpath, ", ") << "') does not correspond to ConsensusXML (run not found). Check input!\n";
          return ILLEGAL_PARAMETERS;
        }
        for (PeptideIdentification& pep_id : new_upep_ids)
        {
          pep_id.setIdentifier(ptr_cmap->second);
        }

        // annotate the RT alignment
        if (qc_rt_alignment.isRunnable(status))
        {
          qc_rt_alignment.compute(new_upep_ids, trafo_descr);
        }


        // save the just calculated IDs for appending to Cmap later (not now, because the vector might resize and invalidate our PepID*).
        all_new_upep_ids.insert(all_new_upep_ids.end(), new_upep_ids.begin(), new_upep_ids.end());
      }

      StringList out_feat = getStringList_("out_feat");
      if (!out_feat.empty())
      {
        FeatureXMLFile().store(out_feat[i], *fmap);
      }
      //------------------------------------------------------------- 
      // Annotate calculated meta values from FeatureMap to given ConsensusMap
      //-------------------------------------------------------------

      // copy MetaValues of unassigned PepIDs
      addPepIDMetaValues_(fmap->getUnassignedPeptideIdentifications(), customID_to_cpepID, mp_f.identifier_to_msrunpath);

      // copy MetaValues of assigned PepIDs
      for (Feature& feature : *fmap)
      {
        addPepIDMetaValues_(feature.getPeptideIdentifications(), customID_to_cpepID, mp_f.identifier_to_msrunpath);
      }
    }
    // mztab writer requires single PIs per CF
    // adds 'feature_id' metavalue to all PIs before moving them to remember the uniqueID of the CF
    IDConflictResolverAlgorithm::resolve(cmap);

    // check if all PepIDs of ConsensusMap appeared in a FeatureMap
    bool incomplete_features {false};
    auto f =
        [&incomplete_features](const PeptideIdentification& pep_id)
        {
          if (!pep_id.getHits().empty() && !pep_id.getHits()[0].metaValueExists("missed_cleavages"))
          {
            OPENMS_LOG_ERROR << "A PeptideIdentification in the ConsensusXML with sequence " << pep_id.getHits()[0].getSequence().toString()
                             << ", RT '" << pep_id.getRT() << "', m/z '" << pep_id.getMZ() << "' and identifier '" << pep_id.getIdentifier()
                             << "' does not appear in any of the given FeatureXMLs. Check your input!\n";
            incomplete_features = true;
          }
        };
    cmap.applyFunctionOnPeptideIDs(f, true);
    if (incomplete_features) return ILLEGAL_PARAMETERS;
    
    // add new PeptideIdentifications (for unidentified MS2 spectra)
    cmap.getUnassignedPeptideIdentifications().insert(cmap.getUnassignedPeptideIdentifications().end(), all_new_upep_ids.begin(), all_new_upep_ids.end());

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    String out_cm = getStringOption_("out_cm");
    if (!out_cm.empty())
    {
      ConsensusXMLFile().store(out_cm, cmap);
    }

    String out = getStringOption_("out");
    if (!out.empty())
    {
      MzTab mztab = MzTab::exportConsensusMapToMzTab(cmap, in_cm, true, true, true, true, "QC export from OpenMS");
      MzTabMetaData meta = mztab.getMetaData();
      qc_tic.addMetaDataMetricsToMzTab(meta);
      qc_ms2ir.addMetaDataMetricsToMzTab(meta);
      mztab.setMetaData(meta);

      MzTabFile mztab_out;
      mztab_out.store(out, mztab);
    }

    return EXECUTION_OK;
  }

private:
  StringList updateFileStatus_(QCBase::Status& status, UInt64& number_exps, const String& port, const QCBase::Requires& req) const
  {
    // since files are optional, leave function if none are provided by the user
    StringList files = getStringList_(port);
    if (!files.empty())
    {
      if (number_exps == 0) number_exps = files.size(); // Number of experiments is determined from first non empty file list.
      if (number_exps != files.size()) // exit if any file list has different length
      {
        throw(Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, port + ": invalid number of files. Expected were " + number_exps + ".\n"));
      }
      status |= req;
    }
    return files;
  }
  
  void fillConsensusPepIDMap_(vector<PeptideIdentification>& cpep_ids,
                              const map<String, StringList>& identifier_to_msrunpath,
                              multimap<String, PeptideIdentification*>& customID_to_cpepID) const
  {
    for (PeptideIdentification& cpep_id : cpep_ids)
    {
      if (!cpep_id.metaValueExists("spectrum_reference"))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Spectrum reference missing at PeptideIdentification.");
      }
      const auto& ms_run_path = identifier_to_msrunpath.at(cpep_id.getIdentifier());

      String UID; //< unique ID to identify the PepID
      if (ms_run_path.size() == 1)
      {
        UID = ms_run_path[0] + cpep_id.getMetaValue("spectrum_reference").toString();
      }
      else if (cpep_id.metaValueExists("map_index"))
      {
        UID = cpep_id.getMetaValue("map_index").toString() + cpep_id.getMetaValue("spectrum_reference").toString();
      }
      else
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Multiple files in a run, but no map_index in PeptideIdentification found.");
      }
      customID_to_cpepID.insert(make_pair(UID, &cpep_id));
    }
  }


  void addPepIDMetaValues_(const vector<PeptideIdentification>& f_pep_ids,
    const multimap<String, PeptideIdentification*>& customID_to_pepID,
    const map<String, StringList>& fidentifier_to_msrunpath) const
  {
    for (const PeptideIdentification& f_pep_id : f_pep_ids)
    {
      // for empty PIs which were created by a metric
      if (f_pep_id.getHits().empty()) continue;

      String UID;
      const auto& ms_run_path = fidentifier_to_msrunpath.at(f_pep_id.getIdentifier());
      if (ms_run_path.size() == 1)
      {
        UID = ms_run_path[0] + f_pep_id.getMetaValue("spectrum_reference").toString();
      }
      else if (f_pep_id.metaValueExists("map_index"))
      {
        UID = f_pep_id.getMetaValue("map_index").toString() + f_pep_id.getMetaValue("spectrum_reference").toString();
      }
      else
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Multiple files in a run, but no map_index in PeptideIdentification found.");
      }

      const auto range = customID_to_pepID.equal_range(UID);

      for (auto it_pep = range.first; it_pep != range.second; ++it_pep) // OMS_CODING_TEST_EXCLUDE
      {
        // copy all MetaValues that are at PepID level
        it_pep->second->addMetaValues(f_pep_id);

        // copy all MetaValues that are at best Hit level
        //TODO check if first = best assumption is met!
        (it_pep->second)->getHits()[0].addMetaValues(f_pep_id.getHits()[0]);
      }
    }
  }
};

// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPQualityControl tool;
  return tool.main(argc, argv);
}

/// @endcond
