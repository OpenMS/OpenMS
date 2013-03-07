// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENTHELPER_H
#define OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENTHELPER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/CVTerm.h>
#include <OpenMS/METADATA/CVTermList.h>

#include <boost/numeric/conversion/cast.hpp>

namespace OpenMS
{
  /**
    @brief This class stores helper structures that are used in multiple
    classes of the TargetedExperiment (e.g. ReactionMonitoringTransition and
    IncludeExcludeTarget).
  */

  namespace TargetedExperimentHelper
  {

    struct Configuration :
      public CVTermList
    {
      String contact_ref;
      String instrument_ref;
      std::vector<CVTermList> validations;

      Configuration & operator=(const Configuration & rhs)
      {
        if (this != &rhs)
        {
          CVTermList::operator=(rhs);
          contact_ref = rhs.contact_ref;
          instrument_ref = rhs.instrument_ref;
          validations = rhs.validations;
        }
        return *this;
      }

    };

    struct CV
    {
      CV(const String & new_id, const String & new_fullname, const String & new_version, const String & new_URI) :
        id(new_id),
        fullname(new_fullname),
        version(new_version),
        URI(new_URI)
      {

      }

      String id;
      String fullname;
      String version;
      String URI;

      bool operator==(const CV & cv) const
      {
        return id == cv.id &&
               fullname == cv.fullname &&
               version == cv.version &&
               URI == cv.URI;
      }

    };

    struct Protein :
      public CVTermList
    {
      Protein() :
        CVTermList()
      {
      }

      String id;
      String sequence;

      bool operator==(const Protein & rhs) const
      {
        return CVTermList::operator==(rhs) &&
               id == rhs.id &&
               sequence == rhs.sequence;
      }

      Protein & operator=(const Protein & rhs)
      {
        if (&rhs != this)
        {
          CVTermList::operator=(rhs);
          id = rhs.id;
          sequence = rhs.sequence;
        }
        return *this;
      }

    };

    class OPENMS_DLLAPI RetentionTime :
      public CVTermList
    {
public:

      RetentionTime() :
        CVTermList()
      {
      }

      RetentionTime(const RetentionTime & rhs) :
        CVTermList(rhs),
        software_ref(rhs.software_ref)
      {
      }

      virtual ~RetentionTime()
      {
      }

      RetentionTime & operator=(const RetentionTime & rhs)
      {
        if (&rhs != this)
        {
          CVTermList::operator=(rhs);
          software_ref = rhs.software_ref;
        }
        return *this;
      }

      bool operator==(const RetentionTime & rhs) const
      {
        return CVTermList::operator==(rhs) &&
               software_ref == rhs.software_ref;
      }

      String software_ref;
    };

    class OPENMS_DLLAPI Compound :
      public CVTermList
    {
public:

      Compound() :
        CVTermList()
      {
      }

      Compound(const Compound & rhs) :
        CVTermList(rhs),
        id(rhs.id),
        rts(rhs.rts)
      {
      }

      Compound & operator=(const Compound & rhs)
      {
        if (this != &rhs)
        {
          CVTermList::operator=(rhs);
          id = rhs.id;
          rts = rhs.rts;
        }
        return *this;
      }

      bool operator==(const Compound & rhs) const
      {
        return CVTermList::operator==(rhs) &&
               id == rhs.id &&
               rts == rhs.rts;
      }

      String id;
      std::vector<RetentionTime> rts;
    };

    class OPENMS_DLLAPI Peptide :
      public CVTermList
    {
public:

      struct Modification :
        public CVTermList
      {
        DoubleReal avg_mass_delta;
        int location;
        DoubleReal mono_mass_delta;
      };

      Peptide() :
        CVTermList()
      {
        charge_ = -1;
        peptide_group_label_ = -1;

        // we store the actual labels in a static vector which allows us to
        // store each label only once.
        static std::vector<String> * init_peptide_group_labels_ = 0;
        if (init_peptide_group_labels_ == 0)
        {
          init_peptide_group_labels_ = new std::vector<String>;
        }
        peptide_group_labels_ = init_peptide_group_labels_;
      }

      Peptide(const Peptide & rhs) :
        CVTermList(rhs),
        rts(rhs.rts),
        id(rhs.id),
        protein_refs(rhs.protein_refs),
        evidence(rhs.evidence),
        sequence(rhs.sequence),
        mods(rhs.mods),
        charge_(rhs.charge_),
        peptide_group_label_(rhs.peptide_group_label_),
        peptide_group_labels_(rhs.peptide_group_labels_)
      {
      }

      Peptide & operator=(const Peptide & rhs)
      {
        if (this != &rhs)
        {
          CVTermList::operator=(rhs);
          rts = rhs.rts;
          id = rhs.id;
          protein_refs = rhs.protein_refs;
          evidence = rhs.evidence;
          sequence = rhs.sequence;
          mods = rhs.mods;
          charge_ = rhs.charge_;
          peptide_group_label_ = rhs.peptide_group_label_;
          peptide_group_labels_ = rhs.peptide_group_labels_;
        }
        return *this;
      }

      bool operator==(const Peptide & rhs) const
      {
        return CVTermList::operator==(rhs) &&
               rts == rhs.rts &&
               id == rhs.id &&
               protein_refs == rhs.protein_refs &&
               evidence == rhs.evidence &&
               sequence == rhs.sequence &&
               mods == rhs.mods &&
               charge_ == rhs.charge_ &&
               peptide_group_label_ == rhs.peptide_group_label_ &&
               peptide_group_labels_ == rhs.peptide_group_labels_;
      }

      void setChargeState(int charge)
      {
        charge_ = charge;
      }

      int getChargeState() const
      {
        return charge_;
      }

      void setPeptideGroupLabel(const String & label)
      {
        for (Size i = 0; i < peptide_group_labels_->size(); i++)
        {
          if ((*peptide_group_labels_)[i] == label)
          {
            peptide_group_label_ = (Int)i;
            return;
          }
        }

        // not found, add it to the list
        peptide_group_label_ = (Int)peptide_group_labels_->size();
        peptide_group_labels_->push_back(label);
      }

      String getPeptideGroupLabel() const
      {
        if (peptide_group_label_ == -1)
        {
          return "";
        }
        return (*peptide_group_labels_)[peptide_group_label_];
      }

      std::vector<RetentionTime> rts;
      String id;
      std::vector<String> protein_refs;
      CVTermList evidence;
      String sequence;
      std::vector<Modification> mods;

protected:
      int charge_;
      Int peptide_group_label_;
      std::vector<String> * peptide_group_labels_;
    };

    struct OPENMS_DLLAPI Contact :
      public CVTermList
    {
      Contact() :
        CVTermList()
      {
      }

      String id;

      bool operator==(const Contact & rhs) const
      {
        return CVTermList::operator==(rhs) &&
               id == rhs.id;
      }

      Contact & operator=(const Contact & rhs)
      {
        if (&rhs != this)
        {
          CVTermList::operator=(rhs);
          id = rhs.id;
        }
        return *this;
      }

    };

    struct OPENMS_DLLAPI Publication :
      public CVTermList
    {
      Publication() :
        CVTermList()
      {
      }

      String id;

      bool operator==(const Publication & rhs) const
      {
        return CVTermList::operator==(rhs) &&
               id == rhs.id;
      }

      Publication & operator=(const Publication & rhs)
      {
        if (&rhs != this)
        {
          CVTermList::operator=(rhs);
          id = rhs.id;
        }
        return *this;
      }

    };

    struct OPENMS_DLLAPI Instrument :
      public CVTermList
    {
      Instrument() :
        CVTermList()
      {
      }

      String id;

      bool operator==(const Instrument & rhs) const
      {
        return CVTermList::operator==(rhs) &&
               id == rhs.id;
      }

      Instrument & operator=(const Instrument & rhs)
      {
        if (&rhs != this)
        {
          CVTermList::operator=(rhs);
          id = rhs.id;
        }
        return *this;
      }

    };

    struct OPENMS_DLLAPI Prediction :
      public CVTermList
    {
      Prediction() :
        CVTermList()
      {
      }

      String software_ref;
      String contact_ref;

      bool operator==(const Prediction & rhs) const
      {
        return CVTermList::operator==(rhs) &&
               contact_ref == rhs.contact_ref &&
               software_ref == rhs.software_ref;
      }

      Prediction & operator=(const Prediction & rhs)
      {
        if (&rhs != this)
        {
          CVTermList::operator=(rhs);
          software_ref = rhs.software_ref;
          contact_ref = rhs.contact_ref;
        }
        return *this;
      }

    };

    struct OPENMS_DLLAPI TraMLProduct :
      public CVTermList
    {
      TraMLProduct() :
        CVTermList()
      {
        charge_ = -1;
      }

      bool operator==(const TraMLProduct & rhs) const
      {
        return CVTermList::operator==(rhs) &&
               charge_ == rhs.charge_ &&
               configuration_list_ == rhs.configuration_list_ &&
               interpretation_list_ == rhs.interpretation_list_;
      }

      TraMLProduct & operator=(const TraMLProduct & rhs)
      {
        if (&rhs != this)
        {
          CVTermList::operator=(rhs);
          charge_ = rhs.charge_;
          configuration_list_ = rhs.configuration_list_;
          interpretation_list_ = rhs.interpretation_list_;
        }
        return *this;
      }

      void setChargeState(int charge)
      {
        charge_ = charge;
      }

      int getChargeState() const
      {
        return charge_;
      }

      const std::vector<Configuration> & getConfigurationList() const
      {
        return configuration_list_;
      }

      void addConfiguration(const Configuration configuration)
      {
        return configuration_list_.push_back(configuration);
      }

      void replaceCVTerms(Map<String, std::vector<CVTerm> > & cv_terms)
      {
        cv_terms_ = cv_terms;
      }

      const std::vector<CVTermList> & getInterpretationList() const
      {
        return interpretation_list_;
      }

      void addInterpretation(const CVTermList interpretation)
      {
        return interpretation_list_.push_back(interpretation);
      }

private:
      int charge_;
      std::vector<Configuration> configuration_list_;
      std::vector<CVTermList> interpretation_list_;

    };

    /// helper function that converts a Peptide object to a AASequence object
    OPENMS_DLLAPI OpenMS::AASequence getAASequence(const Peptide& peptide);

    /// helper function that sets a modifiction on a AASequence object 
    OPENMS_DLLAPI void setModification(int location, int max_size, String modification, OpenMS::AASequence & aas);

  }
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENTHELPER_H
