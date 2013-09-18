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
  @page TOPP_MRMTransitionGroupPicker MRMTransitionGroupPicker

  @brief Picks peaks in MRM chromatograms.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPMRMTransitionGroupPicker 
  : public TOPPBase
{
public:

  TOPPMRMTransitionGroupPicker() 
    : TOPPBase("MRMTransitionGroupPicker", "", false)
  {
  }

protected:

  typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram; // this is the type in which we store the chromatograms for this analysis
  typedef ReactionMonitoringTransition TransitionType;
  typedef TargetedExperiment TargetedExpType;
  typedef MRMTransitionGroup<MSSpectrum <ChromatogramPeak>, TransitionType> MRMTransitionGroupType; // a transition group holds the MSSpectra with the Chromatogram peaks from above

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", StringList::create("mzML"));

    registerInputFile_("tr", "<file>", "", "transition file ('TraML' or 'csv')");
    setValidFormats_("tr", StringList::create("csv,traML"));

    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", StringList::create("featureXML"));

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String &) const
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
      for (AssayMapT::iterator assay_it = assay_map.begin(); assay_it != assay_map.end(); assay_it++)
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
        const TransitionType* transition = assay_map[id][i];
        OpenSwath::ChromatogramPtr cptr = input->getChromatogramById(chromatogram_map[transition->getNativeID()]);
        MSChromatogram<ChromatogramPeak> chromatogram_old;
        OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chromatogram_old, cptr);
        RichPeakChromatogram chromatogram;

        // copy old to new chromatogram
        for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
        {
          ChromatogramPeak peak;
          peak.setMZ(it->getRT());
          peak.setIntensity(it->getIntensity());
          chromatogram.push_back(peak);
        }

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
    FeatureMap<> & output, TargetedExpType& transition_exp)
  {
    MRMTransitionGroupPicker trgroup_picker;
    Param picker_param = getParam_().copy("algorithm:", true);
    trgroup_picker.setParameters(picker_param);

    MRMGroupMapper m;
    m.doMap(input, transition_exp);
    if (!m.allAssaysHaveChromatograms() )
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       "Not all assays could be mapped to chromatograms");
    }

    // Iterating over all the assays
    for (MRMGroupMapper::AssayMapT::iterator assay_it = m.assay_map.begin(); assay_it != m.assay_map.end(); assay_it++)
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
        output.push_back(transition_group.getFeatures()[i]);
      }
    }
  }

  ExitCodes main_(int, const char **)
  {

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String tr_file = getStringOption_("tr");

    boost::shared_ptr<MSExperiment<> > exp ( new MSExperiment<> );
    MzMLFile mzmlfile;
    mzmlfile.setLogType(log_type_);
    mzmlfile.load(in, *exp);

    TargetedExpType transition_exp;
    TraMLFile().load(tr_file, transition_exp);

    FeatureMap<> output;
    OpenSwath::SpectrumAccessPtr input = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
    run_(input, output, transition_exp);

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
