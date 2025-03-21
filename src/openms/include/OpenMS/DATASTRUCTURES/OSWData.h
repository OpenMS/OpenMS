// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>
#include <vector>

namespace OpenMS
{
    class MSExperiment;

    /// hierarchy levels of the OSWData tree
    struct OPENMS_DLLAPI OSWHierarchy
    {
      /// the actual levels
      enum Level
      {
        PROTEIN,
        PEPTIDE,
        FEATURE,
        TRANSITION,
        SIZE_OF_VALUES
      };
      /// strings matching 'Level'
      static const char* LevelName[SIZE_OF_VALUES];
    };

    /// Describes a node in the OSWData model tree. 
    /// If a lower level, e.g. feature, is set, the upper levels need to be set as well.
    /// The lowest level which is set, must be indicated by setting @p lowest.
    struct OPENMS_DLLAPI OSWIndexTrace
    {
      int idx_prot = -1;
      int idx_pep = -1;
      int idx_feat = -1;
      int idx_trans = -1;
      OSWHierarchy::Level lowest = OSWHierarchy::Level::SIZE_OF_VALUES;

      /// is the trace default constructed (=false), or does it point somewhere (=true)?
      bool isSet() const
      {
        return lowest != OSWHierarchy::Level::SIZE_OF_VALUES;
      }

    };

    /// high-level meta data of a transition
    struct OPENMS_DLLAPI OSWTransition
    {
      public:
        /// default c'tor
        OSWTransition() = default;
        /// custom c'tor which fills all the members with data; all members are read-only
        OSWTransition(const String& annotation, const UInt32 id, const float product_mz, const char type, const bool is_decoy);
        OSWTransition(const OSWTransition& rhs) = default;
        OSWTransition& operator=(const OSWTransition& rhs) = default;
        OSWTransition(OSWTransition&& rhs) = default;
        OSWTransition& operator=(OSWTransition&& rhs) = default;
        ~OSWTransition() = default;

        /// e.g. y5/-0.002
        const String& getAnnotation() const
        {
          return annotation_;
        }
        /// ID as used in OSWPeakGroup::transition_ids
        UInt32 getID() const
        {
          return id_;
        }
        /// observed product m/z value
        float getProductMZ() const
        {
          return product_mz_;
        }
        /// b, y
        char getType() const
        {
          return type_;
        }
        /// is this a decoy transition (from a decoy protein/peptide)
        bool isDecoy() const
        {
          return is_decoy_;
        }

      private:
        String annotation_; ///< e.g. y5/-0.002
        UInt32 id_;         ///< ID as used in OSWPeakGroup::transition_ids
        float product_mz_;  ///< observed product m/z value
        char type_;         ///< b, y,
        bool is_decoy_;     ///< is this a decoy transition (from a decoy protein/peptide)
    };

    /**
      A peak group (also called feature) is defined on a small RT range (leftWidth to rightWidth) in a group of extracted transitions (chromatograms).
      The same transitions can be used to defined multiple (usually non-overlapping in RT) peak groups, of which usually only one is correct (lowest q-value).
    */
    class OPENMS_DLLAPI OSWPeakGroup
    {
      public:
        /// fallback value of getQValue() if .osw file did not undergo pyProphet
        static constexpr float QVALUE_MISSING = -1;

        /// just a dummy feature to allow for acceptor output values etc
        OSWPeakGroup() = default;
        
        /// custom c'tor which fills all the members with data; all members are read-only
        OSWPeakGroup(const float rt_experimental, const float rt_left_width, const float rt_right_width, const float rt_delta, std::vector<UInt32>&& transition_ids, const float q_value = -1);
        
        /// Copy c'tor
        OSWPeakGroup(const OSWPeakGroup& rhs) = default;
        
        /// copy assignment
        OSWPeakGroup& operator=(const OSWPeakGroup& rhs) = default;

        /// move c'tor
        OSWPeakGroup(OSWPeakGroup&& rhs) = default;
        
        /// move assignment
        OSWPeakGroup& operator=(OSWPeakGroup&& rhs) = default;

        /// observed RT apex position in seconds of the feature
        float getRTExperimental() const
        {
          return rt_experimental_;
        }
        /// RT position in seconds of the left border
        float getRTLeftWidth() const
        {
          return rt_left_width_;
        }
        /// RT position in seconds of the right border
        float getRTRightWidth() const
        {
          return rt_right_width_;
        }
        /// RT difference in seconds to the expected RT
        float getRTDelta() const
        {
          return rt_delta_;
        }
        /// this might return QVALUE_MISSING if q-value is not annotated in the OSW file
        float getQValue() const
        {
          return q_value_;
        }
        /// get the transition ids (can be mapped to the chromatogram XICs in sqMass data)
        const std::vector<UInt32>& getTransitionIDs() const
        {
          return transition_ids_;
        }

      private:
        float rt_experimental_{ 0 }; ///< rt apex of this feature in seconds (averaged across all transitions)
        float rt_left_width_{ 0 };   ///< rt start in seconds
        float rt_right_width_{ 0 };  ///< rt end in seconds
        float rt_delta_{ 0 };        ///< rt offset from expected distance
        float q_value_{ -1 };        ///< optional Q-value from pyProphet; equals -1 if not set;
        std::vector<UInt32> transition_ids_; /// many features will point to the same transition (but at different RT);
    };

    /**
      @brief A peptide with a charge state

      An OSWProtein has one or more OSWPeptidePrecursors.

      The OSWPeptidePrecursor contains multiple candidate features (peak groups) of type OSWPeakGroup, only one of which is usually true.

    */
    class OPENMS_DLLAPI OSWPeptidePrecursor
    {
      public:
        /// just a dummy feature to allow for acceptor output values etc
        OSWPeptidePrecursor() = default;
        /// custom c'tor which fills all the members with data; all members are read-only
        OSWPeptidePrecursor(const String& seq, const short charge, const bool decoy, const float precursor_mz, std::vector<OSWPeakGroup>&& features);
        /// Copy c'tor
        OSWPeptidePrecursor(const OSWPeptidePrecursor& rhs) = default;
        /// assignment operator
        OSWPeptidePrecursor& operator=(const OSWPeptidePrecursor& rhs) = default;
        /// move c'tor
        OSWPeptidePrecursor(OSWPeptidePrecursor&& rhs) = default;
        /// move assignment operator
        OSWPeptidePrecursor& operator=(OSWPeptidePrecursor&& rhs) = default;

        /// the peptide sequence (incl. mods)
        const String& getSequence() const
        {
          return seq_;
        }
        /// precursor charge
        short getCharge() const
        {
          return charge_;
        }
        /// is this a decoy feature (from a decoy protein)
        bool isDecoy() const
        {
          return decoy_;
        }
        /// m/z of this charged peptide
        float getPCMz() const
        {
          return precursor_mz_;
        }
        /// candidate explanations
        const std::vector<OSWPeakGroup>& getFeatures() const
        {
          return features_;
        }

      private:
        String seq_;
        short charge_{0};
        bool decoy_{false};
        float precursor_mz_{0};
        std::vector<OSWPeakGroup> features_;
    };

    /**
      @brief A Protein is the highest entity and contains one or more peptides which were found/traced.

    */
    class OPENMS_DLLAPI OSWProtein
    {
      public:
        /// just a dummy feature to allow for acceptor output values etc
        OSWProtein() = default;
        /// custom c'tor which fills all the members with data; all members are read-only
        OSWProtein(const String& accession, const Size id, std::vector<OSWPeptidePrecursor>&& peptides);
        /// Copy c'tor
        OSWProtein(const OSWProtein& rhs) = default;
        /// assignment operator
        OSWProtein& operator=(const OSWProtein& rhs) = default;
        /// move c'tor
        OSWProtein(OSWProtein&& rhs) = default;
        /// move assignment operator
        OSWProtein& operator=(OSWProtein&& rhs) = default;

        const String& getAccession() const
        {
          return accession_;
        }

        Size getID() const
        {
          return id_;
        }

        const std::vector<OSWPeptidePrecursor>& getPeptidePrecursors() const
        {
          return peptides_;
        }

      private:
        String accession_;
        Size id_;
        std::vector<OSWPeptidePrecursor> peptides_;
    };


    /**
      @brief Holds all or partial information from an OSW file

      First, fill in all transitions and only then add proteins (which reference transitions via their transition-ids deep down).
      References will be checked and enforced (exception otherwise -- see addProtein()).

    */
    class OPENMS_DLLAPI OSWData
    {
      public:
        /// Adds a transition; do this before adding Proteins
        void addTransition(const OSWTransition& tr)
        {
          transitions_.emplace(tr.getID(), tr);
        }

        void addTransition(OSWTransition&& tr)
        {
          UInt32 id = tr.getID();
          transitions_.emplace(id, std::move(tr));
        }

        /// Adds a protein, which has all its subcomponents already populated
        /// All transition references internally are checked to make sure
        /// they are valid.
        /// You can add stub proteins, by omitting their peptide references.
        /// @throws Exception::Precondition() if transition IDs within protein are unknown
        void addProtein(OSWProtein&& prot);

        /// constant accessor to proteins.
        /// There is no mutable access to prevent accidental violation of invariants (i.e. no matching transitions)
        const std::vector<OSWProtein>& getProteins() const
        {
          return proteins_;
        }

        /// Replace existing protein at position @p index
        /// Note: this is NOT the protein ID, but the index into the internal protein vector. See getProteins()
        /// 
        /// @param index A valid index into the getProteins() vector
        /// @param protein The protein to replace the existing one
        /// @throws Exception::Precondition() if transition IDs within protein are unknown
        void setProtein(const Size index, OSWProtein&& protein)
        {
          checkTransitions_(protein);
          proteins_[index] = std::move(protein);
        }

        /// get the total number of transitions (chromatograms)
        Size transitionCount() const
        {
          return transitions_.size();
        }

        /// obtain a certain transition meta information with @p id (this matches the ID of a chromatogram in an sqMass file).
        const OSWTransition& getTransition(const UInt32 id) const
        {
          return transitions_.at(id);
        }
        
        /// get all transitions mapped by their ID (UInt32)
        const std::map<UInt32, OSWTransition>& getTransitions() const
        {
          return transitions_;
        }

        void setSqlSourceFile(const String& filename)
        {
          source_file_ = filename;
        }

        const String& getSqlSourceFile() const
        {
          return source_file_;
        }

        void setRunID(const UInt64 run_id)
        {
          run_id_ = run_id;
        }

        UInt64 getRunID() const
        {
          return run_id_;
        }

        /// forget all data
        void clear();

        /// only forget protein data
        void clearProteins();

        /**
          @brief Create an internal mapping from the nativeIDs of all chromatograms (extracted by OpenSwathWorkflow (e.g. as sqMass file)) to their index (.getChromatograms[index])

          The mapping is stored internally and can be used to translate transition.ids (which are native_ids) to a chromatogram index of the external sqMass file.

          The mapping can be queried using fromNativeID(int transition.id).

          Make sure that the other OSW data is loaded (at least via OSWFile::readMinimal()) before building this mapping here.

          @param chrom_traces The external sqMass file, which we build the mapping on
          @throws Exception::MissingInformation if any nativeID is not known internally
          @throws Exception::Precondition if the run_ids do not match
        */
        void buildNativeIDResolver(const MSExperiment& chrom_traces);

        /// resolve a transition.id (=nativeID) to a simple chromatogram index (.getChromatograms[index]) of the corresponding sqMass file
        /// Requires prior call to buildNativeIDResolver(), throws Exception::InvalidValue otherwise (or when nativeID is not known)
        UInt fromNativeID(int transition_id) const;

      protected:
        /// All transition references are checked against transitions_ to make sure
        /// they are valid.
        /// @throws Exception::Precondition() if transition IDs within protein are unknown
        void checkTransitions_(const OSWProtein& prot) const;

      private:
        std::map<UInt32, OSWTransition> transitions_;
        std::vector<OSWProtein> proteins_;
        String source_file_;                        ///< remember from which sql OSW file this data is loaded (to lazy load more data)
        UInt64 run_id_;                             ///< the ID of this run from the SQL RUN table
        std::map<UInt32, UInt32> transID_to_index_; ///< map a Transition.ID (==native_id) to a chromatogram index in the sqMass experiment which contains the raw data
    };
    

} // namespace OpenMS