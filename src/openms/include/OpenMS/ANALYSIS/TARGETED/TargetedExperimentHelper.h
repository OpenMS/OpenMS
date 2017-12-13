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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENTHELPER_H
#define OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENTHELPER_H

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/CVTerm.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/METADATA/CVTermListInterface.h>
#include <OpenMS/CHEMISTRY/Residue.h>

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

    /**
      @brief This class stores a retention time structure that is used in TargetedExperiment (representing a TraML file)

      According to the standard, each retention time tag can have one or more
      CV terms describing the retention time in question. The unit and type of
      retention time are stored using the RTUnit and RTType structure while the
      actual value is stored in retention_time_ and can be accessed by getRT /
      setRT. Currently support for RT windows or lower/upper limits is not
      implemented but is available via CV terms.
    */
    class OPENMS_DLLAPI RetentionTime :
      public CVTermListInterface
    {
public:

      enum class RTUnit : std::int8_t
      {
        SECOND = 0,        // RT stored in seconds
        MINUTE,            // RT stored in minutes
        UNKNOWN,           // no stored annotation
        SIZE_OF_RTUNIT
      };

      enum class RTType : std::int8_t
      {
        LOCAL = 0,        // undefined local chromatography
        NORMALIZED,       // standardized reference chromatography
        PREDICTED,        // predicted by referenced software
        HPINS,            // H-PINS "The de facto standard providing the retention times"
        IRT,              // iRT retention time standard
        UNKNOWN,          // no stored annotation
        SIZE_OF_RTTYPE
      };

      RetentionTime() :
        CVTermListInterface(),
        software_ref(""),
        retention_time_unit(RTUnit::SIZE_OF_RTUNIT),
        retention_time_type(RTType::SIZE_OF_RTTYPE),
        retention_time_set_(false),
        retention_time_(0.0)
        // retention_time_width(0.0),
        // retention_time_lower(0.0),
        // retention_time_upper(0.0)
      {
      }

      RetentionTime(const RetentionTime & rhs) :
        CVTermListInterface(rhs),
        software_ref(rhs.software_ref),
        retention_time_unit(rhs.retention_time_unit),
        retention_time_type(rhs.retention_time_type),
        retention_time_set_(rhs.retention_time_set_),
        retention_time_(rhs.retention_time_)
      {
      }

      virtual ~RetentionTime()
      {
      }

      RetentionTime & operator=(const RetentionTime & rhs)
      {
        if (&rhs != this)
        {
          CVTermListInterface::operator=(rhs);
          software_ref = rhs.software_ref;
          retention_time_unit = rhs.retention_time_unit;
          retention_time_type = rhs.retention_time_type;
          retention_time_set_ = rhs.retention_time_set_;
          retention_time_ = rhs.retention_time_;
        }
        return *this;
      }

      bool operator==(const RetentionTime & rhs) const
      {
        return CVTermListInterface::operator==(rhs) &&
               software_ref == rhs.software_ref &&
               retention_time_unit == rhs.retention_time_unit &&
               retention_time_type == rhs.retention_time_type &&
               retention_time_set_ == rhs.retention_time_set_ &&
               retention_time_ == rhs.retention_time_;
      }

      bool isRTset() const 
      {
        return retention_time_set_;
      }
      void setRT(double rt)
      {
        retention_time_ = rt;
        retention_time_set_ = true;
      }
      double getRT() const
      {
        OPENMS_PRECONDITION(isRTset(), "RT needs to be set")
        return retention_time_;
      }

      String software_ref;
      RTUnit retention_time_unit;
      RTType retention_time_type;

private:

      bool retention_time_set_;
      double retention_time_;
      // double retention_time_width;
      // double retention_time_lower;
      // double retention_time_upper;
    };

    class OPENMS_DLLAPI PeptideCompound :
      public CVTermList
    {
public:

      PeptideCompound() :
        CVTermList(),
        charge_(0),
        charge_set_(false)
      {
      }

      PeptideCompound(const PeptideCompound & rhs) :
        CVTermList(rhs),
        id(rhs.id),
        rts(rhs.rts),
        charge_(rhs.charge_),
        charge_set_(rhs.charge_set_)
      {
      }

      PeptideCompound & operator=(const PeptideCompound & rhs)
      {
        if (this != &rhs)
        {
          CVTermList::operator=(rhs);
          rts = rhs.rts;
          id = rhs.id;
          charge_ = rhs.charge_;
          charge_set_ = rhs.charge_set_;
        }
        return *this;
      }

      bool operator==(const PeptideCompound & rhs) const
      {
        return CVTermList::operator==(rhs) &&
               rts == rhs.rts &&
               id == rhs.id &&
               charge_ == rhs.charge_ &&
               charge_set_ == rhs.charge_set_;
      }

      /// Set the peptide or compound charge state
      void setChargeState(int charge)
      {
        charge_ = charge;
        charge_set_ = true;
      }

      /// Whether peptide or compound has set charge state
      bool hasCharge() const
      {
        return charge_set_;
      }

      /// Return the peptide or compound charge state
      int getChargeState() const
      {
        OPENMS_PRECONDITION(charge_set_, "Cannot return charge which was never set")
        return charge_;
      }

      //@{

      /// Get compound or peptide retentiontime
      bool hasRetentionTime() const
      {
        return (!rts.empty() && rts[0].isRTset());
      }

      double getRetentionTime() const
      {
        if (!hasRetentionTime())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "No retention time information available");
        }
        return rts[0].getRT();
      }

      /// Get compound or peptide retentiontime type
      RetentionTime::RTType getRetentionTimeType() const
      {
        if (!hasRetentionTime())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "No retention time information available");
        }
        return rts[0].retention_time_type;
      }

      /// Get compound or peptide retentiontime unit (minute/seconds)
      RetentionTime::RTUnit getRetentionTimeUnit() const
      {
        if (!hasRetentionTime())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "No retention time information available");
        }
        return rts[0].retention_time_unit;
      }
      //@}

      String id;
      std::vector<RetentionTime> rts;

protected:
      int charge_;
      bool charge_set_;

    };

    class OPENMS_DLLAPI Compound :
      public PeptideCompound
    {
public:

      Compound() :
        theoretical_mass(0.0)
      {
      }

      Compound(const Compound & rhs) :
        PeptideCompound(rhs),
        molecular_formula(rhs.molecular_formula),
        smiles_string(rhs.smiles_string),
        theoretical_mass(rhs.theoretical_mass)
      {
      }

      Compound & operator=(const Compound & rhs)
      {
        if (this != &rhs)
        {
          PeptideCompound::operator=(rhs);
          molecular_formula = rhs.molecular_formula;
          smiles_string = rhs.smiles_string;
          theoretical_mass = rhs.theoretical_mass;
        }
        return *this;
      }

      bool operator==(const Compound & rhs) const
      {
        return PeptideCompound::operator==(rhs) &&
               molecular_formula == rhs.molecular_formula &&
               smiles_string == rhs.smiles_string &&
               theoretical_mass == rhs.theoretical_mass;
      }

      String molecular_formula;
      String smiles_string;
      double theoretical_mass;

protected:

    };

    class OPENMS_DLLAPI Peptide :
      public PeptideCompound
    {
public:

      struct Modification :
        public CVTermListInterface
      {
        double avg_mass_delta;
        double mono_mass_delta;
        Int32 location;
        Int32 unimod_id;

        Modification() :
          CVTermListInterface(),
          location(-1),
          unimod_id(-1)
        {
        }

      };

      Peptide() :
        PeptideCompound()
      {
      }

      Peptide(const Peptide & rhs) :
        PeptideCompound(rhs),
        protein_refs(rhs.protein_refs),
        evidence(rhs.evidence),
        sequence(rhs.sequence),
        mods(rhs.mods),
        peptide_group_label_(rhs.peptide_group_label_)
      {
      }

      Peptide & operator=(const Peptide & rhs)
      {
        if (this != &rhs)
        {
          PeptideCompound::operator=(rhs);
          protein_refs = rhs.protein_refs;
          evidence = rhs.evidence;
          sequence = rhs.sequence;
          mods = rhs.mods;
          peptide_group_label_ = rhs.peptide_group_label_;
        }
        return *this;
      }

      bool operator==(const Peptide & rhs) const
      {
        return PeptideCompound::operator==(rhs) &&
               protein_refs == rhs.protein_refs &&
               evidence == rhs.evidence &&
               sequence == rhs.sequence &&
               mods == rhs.mods &&
               peptide_group_label_ == rhs.peptide_group_label_;
      }

      /** @name The peptide group label specifies to non-labeled peptide group to which the peptide belongs
       *
       * MS:1000893: "An arbitrary string label used to mark a set of peptides
       * that belong together in a set, whereby the members are differentiated
       * by different isotopic labels. For example, the heavy and light forms
       * of the same peptide will both be assigned the same peptide group
       * label." [PSI:MS]
       *
     */
      //@{
      /// Set the peptide group label
      void setPeptideGroupLabel(const String & label)
      {
        peptide_group_label_ = label;
      }

      /// Get the peptide group label
      String getPeptideGroupLabel() const
      {
        return peptide_group_label_;
      }
      //@}

      std::vector<String> protein_refs;
      CVTermList evidence;
      String sequence;
      std::vector<Modification> mods;

protected:
      String peptide_group_label_;
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

    struct OPENMS_DLLAPI Interpretation :
      public CVTermListInterface
    {

      /*
      enum ResidueType
      {
        Full = 0,       // with N-terminus and C-terminus
        Internal,       // internal, without any termini
        NTerminal,      // only N-terminus
        CTerminal,      // only C-terminus
        AIon,           // MS:1001229 N-terminus up to the C-alpha/carbonyl carbon bond
        BIon,           // MS:1001224 N-terminus up to the peptide bond
        CIon,           // MS:1001231 N-terminus up to the amide/C-alpha bond
        XIon,           // MS:1001228 amide/C-alpha bond up to the C-terminus
        YIon,           // MS:1001220 peptide bond up to the C-terminus
        ZIon,           // MS:1001230 C-alpha/carbonyl carbon bond
        Precursor,      // MS:1001523 Precursor ion
        BIonMinusH20,   // MS:1001222 b ion without water
        YIonMinusH20,   // MS:1001223 y ion without water
        BIonMinusNH3,   // MS:1001232 b ion without ammonia
        YIonMinusNH3,   // MS:1001233 y ion without ammonia
        NonIdentified,  // MS:1001240 Non-identified ion
        Unannotated,    // no stored annotation
        SizeOfResidueType
      };
      */

      typedef Residue::ResidueType IonType; // Interpretation IonType

      unsigned char ordinal; // MS:1000903 (product ion series ordinal)
      unsigned char rank; // MS:1000926 (product interpretation rank)
      IonType iontype; // which type of ion (b/y/z/ ...), see Residue::ResidueType

      // Constructor
      Interpretation() :
        CVTermListInterface(),
        ordinal(0),
        rank(0),
        iontype(Residue::Unannotated) // Unannotated does not imply any MS OBO term
      {
      }

      // Copy constructor
      Interpretation(const Interpretation & rhs) :
        CVTermListInterface(rhs),
        ordinal(rhs.ordinal),
        rank(rhs.rank),
        iontype(rhs.iontype)
      {
      }

      /** @name Operators assignment, equality, inequality
      */
      //@{
      bool operator==(const Interpretation & rhs) const
      {
        return CVTermListInterface::operator==(rhs) &&
               ordinal == rhs.ordinal &&
               rank == rhs.rank &&
               iontype == rhs.iontype;
      }

      Interpretation & operator=(const Interpretation & rhs)
      {
        if (&rhs != this)
        {
          CVTermListInterface::operator=(rhs);
          ordinal = rhs.ordinal;
          rank = rhs.rank;
          iontype = rhs.iontype;
        }
        return *this;
      }

      bool operator!=(const Interpretation & rhs) const
      {
        return !(operator==(rhs));
      }
      //@}

    };

    struct OPENMS_DLLAPI TraMLProduct :
      public CVTermListInterface
    {
      TraMLProduct() :
        CVTermListInterface(),
        charge_(0),
        charge_set_(false),
        mz_(0)
      {
      }

      bool operator==(const TraMLProduct & rhs) const
      {
        return CVTermListInterface::operator==(rhs) &&
               charge_ == rhs.charge_ &&
               charge_set_ == rhs.charge_set_ &&
               mz_ == rhs.mz_ &&
               configuration_list_ == rhs.configuration_list_ &&
               interpretation_list_ == rhs.interpretation_list_;
      }

      TraMLProduct & operator=(const TraMLProduct & rhs)
      {
        if (&rhs != this)
        {
          CVTermListInterface::operator=(rhs);
          charge_ = rhs.charge_;
          charge_set_ = rhs.charge_set_;
          mz_ = rhs.mz_;
          configuration_list_ = rhs.configuration_list_;
          interpretation_list_ = rhs.interpretation_list_;
        }
        return *this;
      }

      void setChargeState(int charge)
      {
        charge_ = charge;
        charge_set_ = true;
      }

      /// Whether product has set charge state
      bool hasCharge() const
      {
        return charge_set_;
      }

      int getChargeState() const
      {
        OPENMS_PRECONDITION(charge_set_, "Cannot return charge which was never set")
        return charge_;
      }

      double getMZ() const
      {
        return mz_;
      }

      void setMZ(double mz)
      {
        mz_ = mz;
      }

      const std::vector<Configuration> & getConfigurationList() const
      {
        return configuration_list_;
      }

      void addConfiguration(const Configuration configuration)
      {
        return configuration_list_.push_back(configuration);
      }

      const std::vector<Interpretation> & getInterpretationList() const
      {
        return interpretation_list_;
      }

      void addInterpretation(const Interpretation interpretation)
      {
        return interpretation_list_.push_back(interpretation);
      }

      void resetInterpretations()
      {
        return interpretation_list_.clear();
      }

private:
      int charge_;
      bool charge_set_;
      double mz_;
      std::vector<Configuration> configuration_list_;
      std::vector<Interpretation> interpretation_list_;

    };

    /// helper function that converts a Peptide object to a AASequence object
    OPENMS_DLLAPI OpenMS::AASequence getAASequence(const Peptide& peptide);

    /// helper function that sets a modification on a AASequence object
    OPENMS_DLLAPI void setModification(int location, int max_size, String modification, OpenMS::AASequence & aas);

  }

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENTHELPER_H
