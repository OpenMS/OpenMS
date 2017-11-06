
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
//

#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>

using namespace std;

namespace OpenMS
{

  MRMMapping::MRMMapping() :
      DefaultParamHandler("MRMMapping")
    {
      defaults_.setValue("precursor_tolerance", 0.1, "Precursor tolerance when mapping (in Th)");
      defaults_.setValue("product_tolerance", 0.1, "Product tolerance when mapping (in Th)");

      defaults_.setValue("map_multiple_assays", "false", "Allow to map multiple assays to chromatograms and duplicate these chromatograms in the output.");
      defaults_.setValidStrings("map_multiple_assays", ListUtils::create<String>("true,false"));

      defaults_.setValue("error_on_unmapped", "false", "Treat remaining, unmapped chromatograms as an error");
      defaults_.setValidStrings("error_on_unmapped", ListUtils::create<String>("true,false"));

      // write defaults into Param object param_
      defaultsToParam_();
      updateMembers_();
    }

  void MRMMapping::updateMembers_()
  {
    precursor_tol_ = (double)param_.getValue("precursor_tolerance");
    product_tol_ = (double)param_.getValue("product_tolerance");
    map_multiple_assays_ = (bool)param_.getValue("map_multiple_assays").toBool();
    error_on_unmapped_ = (bool)param_.getValue("error_on_unmapped").toBool();
  }

  void MRMMapping::mapExperiment(const OpenMS::PeakMap& chromatogram_map, 
      const OpenMS::TargetedExperiment& targeted_exp,
      OpenMS::PeakMap& output)
  {
    // copy all meta data from old MSExperiment
    output = (ExperimentalSettings)chromatogram_map;
    output.clear(false);
    std::vector<MSChromatogram > empty_chromats;
    output.setChromatograms(empty_chromats);

    int notmapped = 0;
    for (Size i = 0; i < chromatogram_map.getChromatograms().size(); i++)
    {
      // try to find the best matching transition for this chromatogram
      const MSChromatogram& chromatogram = chromatogram_map.getChromatograms()[i];

      bool prec_product_set = !( std::fabs(chromatogram.getPrecursor().getMZ()) < 1e-5 && 
                                 std::fabs(chromatogram.getProduct().getMZ()) < 1e-5);
      if (!prec_product_set)
      {
        if (map_multiple_assays_)
        {
          LOG_WARN << "Warning: Chromatogram " + 
            String(chromatogram.getNativeID()) + " has no precursor or product m/z recorded, mapping may not work." << std::endl;
        }
        else
        {
          LOG_WARN << "Skip mapping for chromatogram " + 
            String(chromatogram.getNativeID()) + " since no precursor or product m/z was recorded." << std::endl;
          continue;
        }
      }

      std::vector<MSChromatogram > mapped_chroms;
      for (Size j = 0; j < targeted_exp.getTransitions().size(); j++)
      {
        if (fabs(chromatogram.getPrecursor().getMZ() - targeted_exp.getTransitions()[j].getPrecursorMZ()) < precursor_tol_ &&
            fabs(chromatogram.getProduct().getMZ()   - targeted_exp.getTransitions()[j].getProductMZ())   < product_tol_)
        {
          LOG_DEBUG << "Mapping chromatogram " << i << " to transition " << j << " (" << targeted_exp.getTransitions()[j].getNativeID() << ")"
             " with precursor mz " << chromatogram.getPrecursor().getMZ() << " / " <<  targeted_exp.getTransitions()[j].getPrecursorMZ() <<
             " and product mz " << chromatogram.getProduct().getMZ() << " / " <<  targeted_exp.getTransitions()[j].getProductMZ() << std::endl;

          // Create precursor and set the peptide sequence
          MSChromatogram c = chromatogram_map.getChromatograms()[i];
          Precursor precursor = c.getPrecursor();
          String pepref = targeted_exp.getTransitions()[j].getPeptideRef();
          for (Size pep_idx = 0; pep_idx < targeted_exp.getPeptides().size(); pep_idx++)
          {
            const OpenMS::TargetedExperiment::Peptide * pep = &targeted_exp.getPeptides()[pep_idx];
            if (pep->id == pepref)
            {
              precursor.setMetaValue("peptide_sequence", pep->sequence);
              break;
            }
          }
          // add precursor to chromatogram
          c.setPrecursor(precursor);

          // Set the id of the chromatogram, using the id of the transition (this gives directly the mapping of the two
          c.setNativeID(targeted_exp.getTransitions()[j].getNativeID());

          mapped_chroms.push_back(c);
        }
      }

      // Check whether we have mapped this chromatogram to at least one transition:
      //  - warn if no mapping occured
      //  - else append all mapped chromatograms (if we allow multiple mappings)
      //  - else append the first mapped chromatograms (if we don't allow multiple mappings)
      if (mapped_chroms.empty())
      {
        LOG_WARN << "Did not find a mapping for chromatogram " + String(i) + " with transition " + String(chromatogram.getPrecursor().getMZ()) + \
          " -> " + String(chromatogram.getProduct().getMZ()) +  "! Maybe try to increase your mapping tolerance." << std::endl;
        notmapped++;
        if (error_on_unmapped_)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Did not find a mapping for chromatogram " + String(i) + "! Maybe try to increase your mapping tolerance.");
        }
      }
      else if (map_multiple_assays_)
      {
        for (auto & c : mapped_chroms) output.addChromatogram(c);
      }
      else
      {
        if (mapped_chroms.size() == 1) output.addChromatogram(mapped_chroms[0]);
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Chromatogram " + String(chromatogram.getNativeID()) + \
           " with " + String(chromatogram.getPrecursor().getMZ()) + \
            " -> " + String(chromatogram.getProduct().getMZ()) + \
              "maps to multiple assays! Either decrease your mapping tolerance or set map_multiple_assays to true.");
        }
      }
    }

    if (notmapped > 0)
    {
      LOG_WARN << "Could not find mapping for " << notmapped  << " chromatogram(s)." << std::endl;
      if (error_on_unmapped_)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Found " + String(notmapped) + \
            " unmapped chromatograms, disable error_on_unmapped to continue.");
      }
    }


  }

} //namespace

