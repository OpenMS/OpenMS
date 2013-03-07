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
// $Maintainer: Hannes Roest $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENT_H
#define OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENT_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h>
#include <OpenMS/METADATA/CVTerm.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief This class stores an prediction of an SRM/MRM transition
  */
  class OPENMS_DLLAPI TargetedExperiment
  {
public:

    typedef TargetedExperimentHelper::CV CV;
    typedef TargetedExperimentHelper::Protein Protein;
    typedef TargetedExperimentHelper::RetentionTime RetentionTime;
    typedef TargetedExperimentHelper::Compound Compound;
    typedef TargetedExperimentHelper::Peptide Peptide;
    typedef TargetedExperimentHelper::Contact Contact;
    typedef TargetedExperimentHelper::Publication Publication;
    typedef TargetedExperimentHelper::Instrument Instrument;
    typedef TargetedExperimentHelper::Prediction Prediction;

    typedef std::map<String, const Protein *> ProteinReferenceMapType;
    typedef std::map<String, const Peptide *> PeptideReferenceMapType;

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    TargetedExperiment();

    /// copy constructor
    TargetedExperiment(const TargetedExperiment & rhs);

    /// destructor
    virtual ~TargetedExperiment();
    //@}

    /// assignment operator
    TargetedExperiment & operator=(const TargetedExperiment & rhs);

    /** @name Predicates
    */
    //@{
    bool operator==(const TargetedExperiment & rhs) const;
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
    TargetedExperiment & operator+=(const TargetedExperiment & rhs);

    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    void clear(bool clear_meta_data);

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

    const std::vector<Protein> & getProteins() const;

    const Protein & getProteinByRef(const String & ref);

    void addProtein(const Protein & protein);

    // compound list
    void setCompounds(const std::vector<Compound> & rhs);

    const std::vector<Compound> & getCompounds() const;

    void addCompound(const Compound & rhs);

    void setPeptides(const std::vector<Peptide> & rhs);

    const std::vector<Peptide> & getPeptides() const;

    const Peptide & getPeptideByRef(const String & ref);

    void addPeptide(const Peptide & rhs);

    /// set transition list
    void setTransitions(const std::vector<ReactionMonitoringTransition> & transitions);

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

protected:

    void createProteinReferenceMap_();

    void createPeptideReferenceMap_();

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

    ProteinReferenceMapType protein_reference_map_;

    bool protein_reference_map_dirty_;

    PeptideReferenceMapType peptide_reference_map_;

    bool peptide_reference_map_dirty_;

  };


  namespace TargetedExperimentHelper
  {
  } // namespace TargetedExperimentHelper


} // namespace OpenMS

#endif // OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENT_H
