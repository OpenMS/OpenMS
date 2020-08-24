// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
  TargetedExperiment::TargetedExperiment() :
    protein_reference_map_dirty_(true),
    peptide_reference_map_dirty_(true),
    compound_reference_map_dirty_(true)
  {
  }

  TargetedExperiment::TargetedExperiment(const TargetedExperiment & rhs) :
    cvs_(rhs.cvs_),
    contacts_(rhs.contacts_),
    publications_(rhs.publications_),
    instruments_(rhs.instruments_),
    targets_(rhs.targets_),
    software_(rhs.software_),
    proteins_(rhs.proteins_),
    compounds_(rhs.compounds_),
    peptides_(rhs.peptides_),
    transitions_(rhs.transitions_),
    include_targets_(rhs.include_targets_),
    exclude_targets_(rhs.exclude_targets_),
    source_files_(rhs.source_files_),
    protein_reference_map_dirty_(true),
    peptide_reference_map_dirty_(true),
    compound_reference_map_dirty_(true)
  {
  }

  TargetedExperiment::~TargetedExperiment()
  {
  }

  TargetedExperiment & TargetedExperiment::operator=(const TargetedExperiment & rhs)
  {
    if (&rhs != this)
    {
      cvs_ = rhs.cvs_;
      contacts_ = rhs.contacts_;
      publications_ = rhs.publications_;
      instruments_ = rhs.instruments_;
      targets_ = rhs.targets_;
      software_ = rhs.software_;
      proteins_ = rhs.proteins_;
      compounds_ = rhs.compounds_;
      peptides_ = rhs.peptides_;
      transitions_ = rhs.transitions_;
      include_targets_ = rhs.include_targets_;
      exclude_targets_ = rhs.exclude_targets_;
      source_files_ = rhs.source_files_;
      protein_reference_map_dirty_ = true;
      peptide_reference_map_dirty_ = true;
      compound_reference_map_dirty_ = true;
    }
    return *this;
  }

  TargetedExperiment TargetedExperiment::operator+(const TargetedExperiment & rhs) const
  {
    TargetedExperiment tmp(*this);
    tmp += rhs;
    return tmp;
  }

  TargetedExperiment & TargetedExperiment::operator+=(const TargetedExperiment & rhs)
  {
    protein_reference_map_dirty_ = true;
    peptide_reference_map_dirty_ = true;
    compound_reference_map_dirty_ = true;

    // merge these:
    cvs_.insert(cvs_.end(), rhs.cvs_.begin(), rhs.cvs_.end());
    contacts_.insert(contacts_.end(), rhs.contacts_.begin(), rhs.contacts_.end());
    publications_.insert(publications_.end(), rhs.publications_.begin(), rhs.publications_.end());
    instruments_.insert(instruments_.end(), rhs.instruments_.begin(), rhs.instruments_.end());
    software_.insert(software_.end(), rhs.software_.begin(), rhs.software_.end());
    proteins_.insert(proteins_.end(), rhs.proteins_.begin(), rhs.proteins_.end());
    compounds_.insert(compounds_.end(), rhs.compounds_.begin(), rhs.compounds_.end());
    peptides_.insert(peptides_.end(), rhs.peptides_.begin(), rhs.peptides_.end());
    transitions_.insert(transitions_.end(), rhs.transitions_.begin(), rhs.transitions_.end());
    include_targets_.insert(include_targets_.end(), rhs.include_targets_.begin(), rhs.include_targets_.end());
    exclude_targets_.insert(exclude_targets_.end(), rhs.exclude_targets_.begin(), rhs.exclude_targets_.end());
    source_files_.insert(source_files_.end(), rhs.source_files_.begin(), rhs.source_files_.end());

    for (Map<String, std::vector<CVTerm> >::const_iterator targ_it = rhs.targets_.getCVTerms().begin(); targ_it != rhs.targets_.getCVTerms().end(); ++targ_it)
    {
      for (std::vector<CVTerm>::const_iterator term_it = targ_it->second.begin(); term_it != targ_it->second.end(); ++term_it)
      {
        targets_.addCVTerm(*term_it);
      }
    }

    // todo: check for double entries
    // transitions, peptides, proteins

    return *this;
  }

  bool TargetedExperiment::operator==(const TargetedExperiment & rhs) const
  {
    return cvs_ == rhs.cvs_ &&
           contacts_ == rhs.contacts_ &&
           publications_ == rhs.publications_ &&
           instruments_ == rhs.instruments_ &&
           targets_ == rhs.targets_ &&
           software_ == rhs.software_ &&
           proteins_ == rhs.proteins_ &&
           compounds_ == rhs.compounds_ &&
           peptides_ == rhs.peptides_ &&
           transitions_ == rhs.transitions_ &&
           include_targets_ == rhs.include_targets_ &&
           exclude_targets_ == rhs.exclude_targets_ &&
           source_files_ == rhs.source_files_;
  }

  bool TargetedExperiment::operator!=(const TargetedExperiment & rhs) const
  {
    return !(operator==(rhs));
  }

  void TargetedExperiment::clear(bool clear_meta_data)
  {
    transitions_.clear();

    if (clear_meta_data)
    {
      cvs_.clear();
      contacts_.clear();
      publications_.clear();
      instruments_.clear();
      targets_ = CVTermList();
      software_.clear();
      proteins_.clear();
      compounds_.clear();
      peptides_.clear();

      include_targets_.clear();
      exclude_targets_.clear();
      source_files_.clear();
      protein_reference_map_.clear();
      peptide_reference_map_.clear();
      compound_reference_map_.clear();

      protein_reference_map_dirty_ = true;
      peptide_reference_map_dirty_ = true;
      compound_reference_map_dirty_ = true;
    }
  }

  void TargetedExperiment::setCVs(const std::vector<CV> & cvs)
  {
    cvs_ = cvs;
  }

  const std::vector<TargetedExperiment::CV> & TargetedExperiment::getCVs() const
  {
    return cvs_;
  }

  void TargetedExperiment::addCV(const CV & cv)
  {
    cvs_.push_back(cv);
  }

  void TargetedExperiment::setContacts(const std::vector<Contact> & contacts)
  {
    contacts_ = contacts;
  }

  const std::vector<TargetedExperiment::Contact> & TargetedExperiment::getContacts() const
  {
    return contacts_;
  }

  void TargetedExperiment::addContact(const Contact & contact)
  {
    contacts_.push_back(contact);
  }

  void TargetedExperiment::setPublications(const std::vector<Publication> & publications)
  {
    publications_ = publications;
  }

  const std::vector<TargetedExperiment::Publication> & TargetedExperiment::getPublications() const
  {
    return publications_;
  }

  void TargetedExperiment::addPublication(const Publication & publication)
  {
    publications_.push_back(publication);
  }

  void TargetedExperiment::setTargetCVTerms(const CVTermList & cv_terms)
  {
    targets_ = cv_terms;
  }

  const CVTermList & TargetedExperiment::getTargetCVTerms() const
  {
    return targets_;
  }

  void TargetedExperiment::addTargetCVTerm(const CVTerm & cv_term)
  {
    targets_.addCVTerm(cv_term);
  }

  void TargetedExperiment::setTargetMetaValue(const String & name, const DataValue & value)
  {
    targets_.setMetaValue(name, value);
  }

  void TargetedExperiment::setInstruments(const std::vector<Instrument> & instruments)
  {
    instruments_ = instruments;
  }

  const std::vector<TargetedExperiment::Instrument> & TargetedExperiment::getInstruments() const
  {
    return instruments_;
  }

  void TargetedExperiment::addInstrument(const Instrument & instrument)
  {
    instruments_.push_back(instrument);
  }

  void TargetedExperiment::setSoftware(const std::vector<Software> & software)
  {
    software_ = software;
  }

  const std::vector<Software> & TargetedExperiment::getSoftware() const
  {
    return software_;
  }

  void TargetedExperiment::addSoftware(const Software & software)
  {
    software_.push_back(software);
  }

  void TargetedExperiment::setProteins(const std::vector<Protein> & proteins)
  {
    protein_reference_map_dirty_ = true;
    proteins_ = proteins;
  }

  const std::vector<TargetedExperiment::Protein> & TargetedExperiment::getProteins() const
  {
    return proteins_;
  }

  const TargetedExperiment::Protein & TargetedExperiment::getProteinByRef(const String & ref) const
  {
    if (protein_reference_map_dirty_)
    {
      createProteinReferenceMap_();
    }
    OPENMS_PRECONDITION(protein_reference_map_.find(ref) != protein_reference_map_.end(), "Could not find protein in map")
    return *(protein_reference_map_[ref]);
  }

  bool TargetedExperiment::hasProtein(const String & ref) const
  {
    if (protein_reference_map_dirty_)
    {
      createProteinReferenceMap_();
    }
    return protein_reference_map_.find(ref) != protein_reference_map_.end();
  }

  void TargetedExperiment::addProtein(const Protein & protein)
  {
    protein_reference_map_dirty_ = true;
    proteins_.push_back(protein);
  }

  void TargetedExperiment::setCompounds(const std::vector<Compound> & compounds)
  {
    compounds_ = compounds;
  }

  const std::vector<TargetedExperiment::Compound> & TargetedExperiment::getCompounds() const
  {
    return compounds_;
  }

  void TargetedExperiment::addCompound(const Compound & rhs)
  {
    compounds_.push_back(rhs);
  }

  void TargetedExperiment::setPeptides(const std::vector<Peptide> & peptides)
  {
    peptide_reference_map_dirty_ = true;
    peptides_ = peptides;
  }

  const std::vector<TargetedExperiment::Peptide> & TargetedExperiment::getPeptides() const
  {
    return peptides_;
  }

  const TargetedExperiment::Peptide & TargetedExperiment::getPeptideByRef(const String & ref) const
  {
    if (peptide_reference_map_dirty_)
    {
      createPeptideReferenceMap_();
    }
    OPENMS_PRECONDITION(hasPeptide(ref), "Cannot return peptide that does not exist, check with hasPeptide() first")
    return *(peptide_reference_map_[ref]);
  }

  const TargetedExperiment::Compound & TargetedExperiment::getCompoundByRef(const String & ref) const
  {
    if (compound_reference_map_dirty_)
    {
      createCompoundReferenceMap_();
    }
    OPENMS_PRECONDITION(hasCompound(ref), "Cannot return compound that does not exist, check with hasCompound() first")
    return *(compound_reference_map_[ref]);
  }

  bool TargetedExperiment::hasPeptide(const String & ref) const
  {
    if (peptide_reference_map_dirty_)
    {
      createPeptideReferenceMap_();
    }
    return peptide_reference_map_.find(ref) != peptide_reference_map_.end();
  }

  bool TargetedExperiment::hasCompound(const String & ref) const
  {
    if (compound_reference_map_dirty_)
    {
      createCompoundReferenceMap_();
    }
    return compound_reference_map_.find(ref) != compound_reference_map_.end();
  }

  void TargetedExperiment::addPeptide(const Peptide & rhs)
  {
    peptide_reference_map_dirty_ = true;
    peptides_.push_back(rhs);
  }

  void TargetedExperiment::setTransitions(const std::vector<ReactionMonitoringTransition> & transitions)
  {
    transitions_ = transitions;
  }

  const std::vector<ReactionMonitoringTransition> & TargetedExperiment::getTransitions() const
  {
    return transitions_;
  }

  void TargetedExperiment::addTransition(const ReactionMonitoringTransition & transition)
  {
    transitions_.push_back(transition);
  }

  void TargetedExperiment::setIncludeTargets(const std::vector<IncludeExcludeTarget> & targets)
  {
    include_targets_ = targets;
  }

  const std::vector<IncludeExcludeTarget> & TargetedExperiment::getIncludeTargets() const
  {
    return include_targets_;
  }

  void TargetedExperiment::addIncludeTarget(const IncludeExcludeTarget & target)
  {
    include_targets_.push_back(target);
  }

  void TargetedExperiment::setExcludeTargets(const std::vector<IncludeExcludeTarget> & targets)
  {
    exclude_targets_ = targets;
  }

  const std::vector<IncludeExcludeTarget> & TargetedExperiment::getExcludeTargets() const
  {
    return exclude_targets_;
  }

  void TargetedExperiment::addExcludeTarget(const IncludeExcludeTarget & target)
  {
    exclude_targets_.push_back(target);
  }

  void TargetedExperiment::setSourceFiles(const std::vector<SourceFile> & source_files)
  {
    source_files_ = source_files;
  }

  const std::vector<SourceFile> & TargetedExperiment::getSourceFiles() const
  {
    return source_files_;
  }

  void TargetedExperiment::addSourceFile(const SourceFile & source_file)
  {
    source_files_.push_back(source_file);
  }

  void TargetedExperiment::sortTransitionsByProductMZ()
  {
    std::sort(transitions_.begin(), transitions_.end(), ReactionMonitoringTransition::ProductMZLess());
  }

  bool TargetedExperiment::containsInvalidReferences() const
  {
    typedef std::vector<OpenMS::TargetedExperiment::Protein> ProteinVectorType;
    typedef std::vector<OpenMS::TargetedExperiment::Peptide> PeptideVectorType;
    typedef std::vector<OpenMS::TargetedExperiment::Compound> CompoundVectorType;
    typedef std::vector<OpenMS::ReactionMonitoringTransition> TransitionVectorType;

    // check that all proteins ids are unique
    std::map<String, int> unique_protein_map;
    for (ProteinVectorType::const_iterator prot_it = getProteins().begin(); prot_it != getProteins().end(); ++prot_it)
    {
      // Create new transition group if it does not yet exist
      if (unique_protein_map.find(prot_it->id) != unique_protein_map.end())
      {
        OPENMS_LOG_ERROR << "Found duplicate protein id (must be unique): " + String(prot_it->id) << std::endl;
        return true;
      }
      unique_protein_map[prot_it->id] = 0;
    }

    // check that all peptide ids are unique
    std::map<String, int> unique_peptide_map;
    for (PeptideVectorType::const_iterator pep_it = getPeptides().begin(); pep_it != getPeptides().end(); ++pep_it)
    {
      // Create new transition group if it does not yet exist
      if (unique_peptide_map.find(pep_it->id) != unique_peptide_map.end())
      {
        OPENMS_LOG_ERROR << "Found duplicate peptide id (must be unique): " + String(pep_it->id) << std::endl;
        return true;
      }
      unique_peptide_map[pep_it->id] = 0;
    }

    // check that all compound ids are unique
    std::map<String, int> unique_compounds_map;
    for (CompoundVectorType::const_iterator comp_it = getCompounds().begin(); comp_it != getCompounds().end(); ++comp_it)
    {
      // Create new transition group if it does not yet exist
      if (unique_compounds_map.find(comp_it->id) != unique_compounds_map.end())
      {
        OPENMS_LOG_ERROR << "Found duplicate compound id (must be unique): " + String(comp_it->id) << std::endl;
        return true;
      }
      unique_compounds_map[comp_it->id] = 0;
    }

    // check that all transition ids are unique
    std::map<String, int> unique_transition_map;
    for (TransitionVectorType::const_iterator tr_it = getTransitions().begin(); tr_it != getTransitions().end(); ++tr_it)
    {
      // Create new transition group if it does not yet exist
      if (unique_transition_map.find(tr_it->getNativeID()) != unique_transition_map.end())
      {
        OPENMS_LOG_ERROR << "Found duplicate transition id (must be unique): " + String(tr_it->getNativeID()) << std::endl;
        return true;
      }
      unique_transition_map[tr_it->getNativeID()] = 0;
    }

    // Check that each peptide has only valid proteins
    for (Size i = 0; i < getPeptides().size(); i++)
    {
      for (std::vector<String>::const_iterator prot_it = getPeptides()[i].protein_refs.begin(); prot_it != getPeptides()[i].protein_refs.end(); ++prot_it)
      {
        if (unique_protein_map.find(*prot_it) == unique_protein_map.end()) 
        {
          OPENMS_LOG_ERROR << "Protein " << *prot_it << " is not present in the provided data structure." << std::endl;
          return true;
        }
      }
    }

    // Check that each peptide has only valid references to peptides and compounds
    for (Size i = 0; i < getTransitions().size(); i++)
    {
      const ReactionMonitoringTransition& tr = getTransitions()[i];
      if (!tr.getPeptideRef().empty())
      {
        if (unique_peptide_map.find(tr.getPeptideRef()) == unique_peptide_map.end()) 
        {
          OPENMS_LOG_ERROR << "Peptide " << tr.getPeptideRef() << " is not present in the provided data structure." << std::endl;
          return true;
        }
      }
      else if (!tr.getCompoundRef().empty())
      {
        if (unique_compounds_map.find(tr.getCompoundRef()) == unique_compounds_map.end()) 
        {
          OPENMS_LOG_ERROR << "Compound " << tr.getPeptideRef() << " is not present in the provided data structure." << std::endl;
          return true;
        }
      }
      else
      {
        // It seems that having no associated compound or peptide is valid as both attributes are optional.
        OPENMS_LOG_WARN << "Transition " << tr.getNativeID() << " does not have a compound or peptide associated with it." << std::endl;
        // return true;
      }
    }
    return false;
  }

  void TargetedExperiment::createProteinReferenceMap_() const
  {
    for (Size i = 0; i < getProteins().size(); i++)
    {
      protein_reference_map_[getProteins()[i].id] = &getProteins()[i];
    }
    protein_reference_map_dirty_ = false;
  }

  void TargetedExperiment::createPeptideReferenceMap_() const
  {
    for (Size i = 0; i < getPeptides().size(); i++)
    {
      peptide_reference_map_[getPeptides()[i].id] = &getPeptides()[i];
    }
    peptide_reference_map_dirty_ = false;
  }

  void TargetedExperiment::createCompoundReferenceMap_() const
  {
    for (Size i = 0; i < getCompounds().size(); i++)
    {
      compound_reference_map_[getCompounds()[i].id] = &getCompounds()[i];
    }
    compound_reference_map_dirty_ = false;
  }

} // namespace OpenMS
