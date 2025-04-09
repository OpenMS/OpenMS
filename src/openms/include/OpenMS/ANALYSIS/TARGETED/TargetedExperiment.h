// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h>
#include <OpenMS/METADATA/CVTerm.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief A description of a targeted experiment containing precursor and production ions.

    A targeted experiment contains transitions used in SRM/MRM as well as
    SWATH-MS/DIA analysis using a targeted approach. This container holds
    descriptions of the precursors and product ions analyzed in such a targeted
    experiment. Generally, the precursor ions can be peptides or small
    molecules (for metabolomics) and each precursor has a set of product ions
    associated with it.

    The TargetedExperiment can be stored to disk either in .traml format using
    the @ref TraMLFile "TraMLFile" or in .tsv format using the TransitionTSVFile.

  */
  class OPENMS_DLLAPI TargetedExperiment
  {
public:
    
    struct OPENMS_DLLAPI SummaryStatistics
    {
      Size protein_count;
      Size peptide_count;
      Size compound_count;
      Size transition_count;
      std::map<ReactionMonitoringTransition::DecoyTransitionType, size_t> decoy_counts; ///< # target/decoy transitions
      bool contains_invalid_references;
    };


    typedef TargetedExperimentHelper::CV CV;
    typedef TargetedExperimentHelper::Protein Protein;
    typedef TargetedExperimentHelper::RetentionTime RetentionTime;
    typedef TargetedExperimentHelper::Compound Compound;
    typedef TargetedExperimentHelper::Peptide Peptide;
    typedef TargetedExperimentHelper::Contact Contact;
    typedef TargetedExperimentHelper::Publication Publication;
    typedef TargetedExperimentHelper::Instrument Instrument;
    typedef TargetedExperimentHelper::Prediction Prediction;
    typedef TargetedExperimentHelper::Interpretation Interpretation;
    typedef ReactionMonitoringTransition Transition;
    typedef Residue IonType; // IonType enum of Interpretation class

    typedef std::map<String, const Protein *> ProteinReferenceMapType;
    typedef std::map<String, const Peptide *> PeptideReferenceMapType;
    typedef std::map<String, const Compound *> CompoundReferenceMapType;

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    TargetedExperiment();

    /// copy constructor
    TargetedExperiment(const TargetedExperiment & rhs);

    /// move constructor
    TargetedExperiment(TargetedExperiment && rhs) noexcept;

    /// destructor
    virtual ~TargetedExperiment();
    //@}

    /// assignment operator
    TargetedExperiment & operator=(const TargetedExperiment & rhs);

    /// move assignment operator
    TargetedExperiment & operator=(TargetedExperiment && rhs) noexcept;

    /** @name Predicates
    */
    //@{
    bool operator==(const TargetedExperiment & rhs) const;

    bool operator!=(const TargetedExperiment & rhs) const;
    //@}

    /**
      @brief Joins two targeted experiments.

      Proteins, peptides and transitions are merged (see operator+= for details).
    */
    TargetedExperiment operator+(const TargetedExperiment & rhs) const;

    /**
      @brief Add one targeted experiment to another.

      @param rhs The targeted experiment to add to this one.
    */
    TargetedExperiment& operator+=(const TargetedExperiment & rhs);
    TargetedExperiment& operator+=(TargetedExperiment && rhs);

    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    void clear(bool clear_meta_data);

    /// return summary stats about this TE.
    SummaryStatistics getSummary() const;

    /** @name Accessors
    */
    //@{
    // cv list
    void setCVs(const std::vector<CV> & cvs);

    const std::vector<CV> & getCVs() const;

    void addCV(const CV & cv);

    // contact list
    void setContacts(const std::vector<Contact> & contacts);

    const std::vector<Contact> & getContacts() const;

    void addContact(const Contact & contact);

    // publication list
    void setPublications(const std::vector<Publication> & publications);

    const std::vector<Publication> & getPublications() const;

    void addPublication(const Publication & publication);

    // target list
    void setTargetCVTerms(const CVTermList & cv_terms);

    const CVTermList & getTargetCVTerms() const;

    void addTargetCVTerm(const CVTerm & cv_term);

    void setTargetMetaValue(const String & name, const DataValue & value);

    // instrument list
    void setInstruments(const std::vector<Instrument> & instruments);

    const std::vector<Instrument> & getInstruments() const;

    void addInstrument(const Instrument & instrument);

    // software list
    void setSoftware(const std::vector<Software> & software);

    const std::vector<Software> & getSoftware() const;

    void addSoftware(const Software & software);

    // protein list
    void setProteins(const std::vector<Protein> & proteins);
    void setProteins(std::vector<Protein> && proteins);

    const std::vector<Protein> & getProteins() const;

    const Protein & getProteinByRef(const String & ref) const;

    bool hasProtein(const String & ref) const;

    void addProtein(const Protein & protein);

    // compound list
    void setCompounds(const std::vector<Compound> & rhs);

    const std::vector<Compound> & getCompounds() const;

    void addCompound(const Compound & rhs);

    void setPeptides(const std::vector<Peptide> & rhs);
    void setPeptides(std::vector<Peptide> && rhs);

    const std::vector<Peptide> & getPeptides() const;

    bool hasPeptide(const String & ref) const;

    const Peptide & getPeptideByRef(const String & ref) const;

    bool hasCompound(const String & ref) const;

    const Compound & getCompoundByRef(const String & ref) const;

    void addPeptide(const Peptide & rhs);

    /// set transition list
    void setTransitions(const std::vector<ReactionMonitoringTransition> & transitions);
    void setTransitions(std::vector<ReactionMonitoringTransition> && transitions);

    /// returns the transition list
    const std::vector<ReactionMonitoringTransition> & getTransitions() const;

    /// adds a transition to the list
    void addTransition(const ReactionMonitoringTransition & transition);

    void setIncludeTargets(const std::vector<IncludeExcludeTarget> & targets);

    const std::vector<IncludeExcludeTarget> & getIncludeTargets() const;

    void addIncludeTarget(const IncludeExcludeTarget & target);

    void setExcludeTargets(const std::vector<IncludeExcludeTarget> & targets);

    const std::vector<IncludeExcludeTarget> & getExcludeTargets() const;

    void addExcludeTarget(const IncludeExcludeTarget & target);

    /// sets the source files
    void setSourceFiles(const std::vector<SourceFile> & source_files);

    /// returns the source file list
    const std::vector<SourceFile> & getSourceFiles() const;

    /// adds a source file to the list
    void addSourceFile(const SourceFile & source_file);
    //@}

    ///@name Sorting peaks
    //@{
    /**
      @brief Lexicographically sorts the transitions by their product m/z.
    */
    void sortTransitionsByProductMZ();
    //@}

    ///@name Sorting peaks
    //@{
    /**
      @brief Lexicographically sorts the transitions by their name.
    */
    void sortTransitionsByName();
    //@}

    /**
      @brief Checks whether the data structure (and the underlying TraML file) contains invalid references

      First checks whether all of the references are unique (protein, peptide,
      compound). Secondly, checks that each reference is valid and points
      either to a protein, peptide or compound. 

      Returns false if the file is valid.
    */
    bool containsInvalidReferences() const;

protected:

    void createProteinReferenceMap_() const;

    void createPeptideReferenceMap_() const;

    void createCompoundReferenceMap_() const;

    std::vector<CV> cvs_;

    std::vector<Contact> contacts_;

    std::vector<Publication> publications_;

    std::vector<Instrument> instruments_;

    CVTermList targets_;

    std::vector<Software> software_;

    std::vector<Protein> proteins_;

    std::vector<Compound> compounds_;

    std::vector<Peptide> peptides_;

    std::vector<ReactionMonitoringTransition> transitions_;

    std::vector<IncludeExcludeTarget> include_targets_;

    std::vector<IncludeExcludeTarget> exclude_targets_;

    std::vector<SourceFile> source_files_;

    mutable ProteinReferenceMapType protein_reference_map_;

    mutable bool protein_reference_map_dirty_;

    mutable PeptideReferenceMapType peptide_reference_map_;

    mutable bool peptide_reference_map_dirty_;

    mutable CompoundReferenceMapType compound_reference_map_;

    mutable bool compound_reference_map_dirty_;

  };

  namespace TargetedExperimentHelper
  {
  } // namespace TargetedExperimentHelper

  /// prints out the summary statistics
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const TargetedExperiment::SummaryStatistics& s);


} // namespace OpenMS

