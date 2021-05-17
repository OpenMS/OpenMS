// --------------------------------------------------------------------------
//           OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the
//  documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/OSWData.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
  const char* OSWHierarchy::LevelName[] = { "protein", "peptide", "feature/peakgroup", "transition" };

  void OSWData::addProtein(OSWProtein&& prot)
  {
    // check if transitions are known
    checkTransitions_(prot);
    proteins_.push_back(std::move(prot));
  }

  void OSWData::clear()
  {
    transitions_.clear();
    proteins_.clear();
  }

  void OSWData::clearProteins()
  {
    proteins_.clear();
  }

  void OSWData::buildNativeIDResolver(const MSExperiment& chrom_traces)
  {
    // first check if the MSExperiment originates from the same run by checking for matching run-ids
    if (chrom_traces.getSqlRunID() != getRunID())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                    "The RUN.ID of the sqMass/MSExperiment ('" + String(chrom_traces.getSqlRunID()) + 
                                    "') and the OSW file ('" + String(getRunID()) + "') does not match. "
                                    "Please use a recent version of OpenSwathWorkflow to create matching data.");
    }
    
    Size chrom_count = chrom_traces.getChromatograms().size();
    for (Size i = 0; i < chrom_count; ++i)
    {
      const auto& chrom = chrom_traces.getChromatograms()[i];
      UInt32 nid;
      try
      {
        nid = chrom.getNativeID().toInt();
      }
      catch (...)
      {
        // probably a precursor native ID, e.g. 5543_precursor_i0 .. currently not handled.
        continue;
      }
      if (transitions_.find(nid) == transitions_.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Transition with nativeID " + (String(nid)) + " not found in OSW data. Make sure the OSW data was loaded!");
      }
      transID_to_index_[nid] = (UInt32)i;
    }
  }

  UInt OSWData::fromNativeID(int transition_id) const
  {
    auto it = transID_to_index_.find(transition_id);
    if (it == transID_to_index_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Native ID not found in sqMass file. Did you load the correct file (corresponding sqMass + OSW file)?", String(transition_id));
    }
    return it->second;
  }

  void OSWData::checkTransitions_(const OSWProtein& prot) const
  {
    for (const auto& pc : prot.getPeptidePrecursors())
    {
      for (const auto& f : pc.getFeatures())
      {
        for (const auto& tr : f.getTransitionIDs())
        {
          if (transitions_.find(tr) == transitions_.end())
          {
            throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Transition with ID " + String(tr) + " was referenced in Protein/Precursor/Feature but is not known!");
          }
        }
      }
    }
  }

  OSWProtein::OSWProtein(const String& accession, const Size id, std::vector<OSWPeptidePrecursor>&& peptides)
    : accession_(accession),
    id_(id),
    peptides_(std::move(peptides))
  {}

  OSWPeptidePrecursor::OSWPeptidePrecursor(const String& seq, const short charge, const bool decoy, const float precursor_mz, std::vector<OSWPeakGroup>&& features)
    : seq_(seq),
    charge_(charge),
    decoy_(decoy),
    precursor_mz_(precursor_mz),
    features_(std::move(features))
  {
  }

  OSWPeakGroup::OSWPeakGroup(const float rt_experimental, const float rt_left_width, const float rt_right_width, const float rt_delta, std::vector<UInt32>&& transition_ids, const float q_value)
    : rt_experimental_(rt_experimental),
    rt_left_width_(rt_left_width),
    rt_right_width_(rt_right_width),
    rt_delta_(rt_delta),
    q_value_(q_value),
    transition_ids_(std::move(transition_ids))
  {
  }

  OSWTransition::OSWTransition(const String& annotation, const UInt32 id, const float product_mz, const char type, const bool is_decoy)
    : annotation_(annotation),
    id_(id),
    product_mz_(product_mz),
    type_(type),
    is_decoy_(is_decoy)
  {}
}
