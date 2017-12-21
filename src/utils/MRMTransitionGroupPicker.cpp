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
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

// files
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

// interfaces
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

// helpers
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>

#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>

using namespace std;

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_MRMTransitionGroupPicker MRMTransitionGroupPicker

  @brief Picks peaks in SRM/MRM chromatograms that belong to the same precursors.

    <CENTER>
        <table>
            <tr>
                <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
                <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ MRMTransitionGroupPicker \f$ \longrightarrow \f$</td>
                <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
            </tr>
            <tr>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathChromatogramExtractor </td>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_OpenSwathFeatureXMLToTSV </td>
            </tr>
            <tr>
                <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MRMMapper </td>
            </tr>
        </table>
    </CENTER>


  This tools accepts a set of chromatograms and picks peaks in them, correctly
  grouping related transitions from the same precursor together. It will
  perform the following steps:
  - Step 1: find features (peaks) in individual chromatograms </li>
  - Step 2: merge these features to consensus features that span multiple chromatograms </li>

  Step 1 is performed by smoothing the individual chromatogram and applying the
  PeakPickerHiRes.

  Step 2 is performed by finding the largest peak overall and use this to
  create a feature, propagating this through all chromatograms.

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_MRMTransitionGroupPicker.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_MRMTransitionGroupPicker.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPMRMTransitionGroupPicker 
  : public TOPPBase
{
public:

  TOPPMRMTransitionGroupPicker() 
    : TOPPBase("MRMTransitionGroupPicker", "Picks peaks in SRM/MRM chromatograms.", false)
  {
  }

protected:

  typedef ReactionMonitoringTransition TransitionType;
  typedef TargetedExperiment TargetedExpType;
  typedef MRMTransitionGroup<MSChromatogram, TransitionType> MRMTransitionGroupType;

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("tr", "<file>", "", "transition file ('TraML' or 'csv')");
    setValidFormats_("tr", ListUtils::create<String>("csv,traML"));

    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String &) const override
  {
    return MRMTransitionGroupPicker().getDefaults();
  }

  struct MRMGroupMapper 
  {
    typedef std::map<String, std::vector< const TransitionType* > > AssayMapT;

    // chromatogram map
    std::map<String, int> chromatogram_map;
    // Map peptide id
    std::map<String, int> assay_peptide_map;
    // Group transitions
    AssayMapT assay_map;

    /// Create the mapping
    void doMap(OpenSwath::SpectrumAccessPtr input, TargetedExpType& transition_exp)
    {
      for (Size i = 0; i < input->getNrChromatograms(); i++)
      {
        chromatogram_map[input->getChromatogramNativeID(i)] = boost::numeric_cast<int>(i);
      }
      for (Size i = 0; i < transition_exp.getPeptides().size(); i++)
      {
        assay_peptide_map[transition_exp.getPeptides()[i].id] = boost::numeric_cast<int>(i);
      }
      for (Size i = 0; i < transition_exp.getTransitions().size(); i++)
      {
        assay_map[transition_exp.getTransitions()[i].getPeptideRef()].push_back(&transition_exp.getTransitions()[i]);
      }
    }

    /// Check that all assays have a corresponding chromatogram
    bool allAssaysHaveChromatograms()
    {
      for (AssayMapT::iterator assay_it = assay_map.begin(); assay_it != assay_map.end(); ++assay_it)
      {
        for (Size i = 0; i < assay_it->second.size(); i++)
        {
          if (chromatogram_map.find(assay_it->second[i]->getNativeID()) == chromatogram_map.end())
          {
            return false;
          }
        }
      }
      return true;
    }

    /// Fill up transition group with paired Transitions and Chromatograms
    void getTransitionGroup(OpenSwath::SpectrumAccessPtr input, MRMTransitionGroupType& transition_group, String id)
    {
      transition_group.setTransitionGroupID(id);

      // Go through all transitions
      for (Size i = 0; i < assay_map[id].size(); i++)
      {

        // Check first whether we have a mapping (e.g. see -force option)
        const TransitionType* transition = assay_map[id][i];
        if (chromatogram_map.find(transition->getNativeID()) == chromatogram_map.end())
        {
          LOG_DEBUG << "Found no matching chromatogram for id " << transition->getNativeID() << std::endl;
          continue;
        }

        OpenSwath::ChromatogramPtr cptr = input->getChromatogramById(chromatogram_map[transition->getNativeID()]);
        MSChromatogram chromatogram;
        OpenSwathDataAccessHelper::convertToOpenMSChromatogram(cptr, chromatogram);

        chromatogram.setMetaValue("product_mz", transition->getProductMZ());
        chromatogram.setMetaValue("precursor_mz", transition->getPrecursorMZ());
        chromatogram.setNativeID(transition->getNativeID());

        // Now add the transition and the chromatogram to the group
        transition_group.addTransition(*transition, transition->getNativeID());
        transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());
      }
    }
    
  };

  void run_(OpenSwath::SpectrumAccessPtr input,
    FeatureMap & output, TargetedExpType& transition_exp, bool force)
  {
    MRMTransitionGroupPicker trgroup_picker;
    Param picker_param = getParam_().copy("algorithm:", true);
    trgroup_picker.setParameters(picker_param);

    MRMGroupMapper m;
    m.doMap(input, transition_exp);
    if (!m.allAssaysHaveChromatograms() && !force)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Not all assays could be mapped to chromatograms");
    }

    // Iterating over all the assays
    for (MRMGroupMapper::AssayMapT::iterator assay_it = m.assay_map.begin(); assay_it != m.assay_map.end(); ++assay_it)
    {
      String id = assay_it->first;

      // Create new transition group if there is none for this peptide
      MRMTransitionGroupType transition_group;
      m.getTransitionGroup(input, transition_group, id);

      // Process the transition_group
      trgroup_picker.pickTransitionGroup(transition_group);

      // Add to output
      for (Size i = 0; i < transition_group.getFeatures().size(); i++)
      {
        MRMFeature mrmfeature = transition_group.getFeatures()[i];
        // Prepare the subordinates for the mrmfeature (process all current
        // features and then append all precursor subordinate features)
        std::vector<Feature> allFeatures = mrmfeature.getFeatures();
        for (std::vector<Feature>::iterator f_it = allFeatures.begin(); f_it != allFeatures.end(); ++f_it)
        {
          f_it->getConvexHulls().clear();
          f_it->ensureUniqueId();
        }
        mrmfeature.setSubordinates(allFeatures); // add all the subfeatures as subordinates
        output.push_back(mrmfeature);
      }
    }
  }

  ExitCodes main_(int, const char **) override
  {

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String tr_file = getStringOption_("tr");
    bool force = getFlag_("force");

    boost::shared_ptr<PeakMap > exp ( new PeakMap );
    MzMLFile mzmlfile;
    mzmlfile.setLogType(log_type_);
    mzmlfile.load(in, *exp);

    TargetedExpType transition_exp;
    TraMLFile().load(tr_file, transition_exp);

    FeatureMap output;
    OpenSwath::SpectrumAccessPtr input = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
    run_(input, output, transition_exp, force);

    output.ensureUniqueId();
    StringList ms_runs;
    exp->getPrimaryMSRunPath(ms_runs);
    output.setPrimaryMSRunPath(ms_runs);
    FeatureXMLFile().store(out, output);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPMRMTransitionGroupPicker tool;
  return tool.main(argc, argv);
}

/// @endcond
