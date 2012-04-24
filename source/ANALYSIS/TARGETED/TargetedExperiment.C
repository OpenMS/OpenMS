// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

#include <algorithm>

namespace OpenMS
{
  TargetedExperiment::TargetedExperiment()
  {
  }

  TargetedExperiment::TargetedExperiment(const TargetedExperiment& rhs)
    : cvs_(rhs.cvs_),
      contacts_(rhs.contacts_),
      publications_(rhs.publications_),
      instruments_(rhs.instruments_),
      software_(rhs.software_),
      proteins_(rhs.proteins_),
      compounds_(rhs.compounds_),
      peptides_(rhs.peptides_),
      transitions_(rhs.transitions_),
      protein_reference_map_dirty_(true),
      peptide_reference_map_dirty_(true)
  {
  }

  TargetedExperiment::~TargetedExperiment()
  {
  }

  TargetedExperiment& TargetedExperiment::operator = (const TargetedExperiment& rhs)
  {
    if (&rhs != this)
    {
      cvs_ = rhs.cvs_;
      contacts_ = rhs.contacts_;
      publications_ = rhs.publications_;
      instruments_ = rhs.instruments_;
      software_ = rhs.software_;
      proteins_ = rhs.proteins_;
      compounds_ = rhs.compounds_;
      peptides_ = rhs.peptides_;
      transitions_ = rhs.transitions_;
      protein_reference_map_dirty_ = true;
      peptide_reference_map_dirty_ = true;
    }
    return *this;
  }

  bool TargetedExperiment::operator == (const TargetedExperiment& rhs) const
  {
    return   cvs_ == rhs.cvs_ &&
            contacts_ == rhs.contacts_ &&
            publications_ == rhs.publications_ &&
            instruments_ == rhs.instruments_ &&
            software_ == rhs.software_ &&
            proteins_ == rhs.proteins_ &&
            compounds_ == rhs.compounds_ &&
            peptides_ == rhs.peptides_ &&
            transitions_ == rhs.transitions_;
  }

  void TargetedExperiment::setCVs(const std::vector<CV>& cvs)
  {
    cvs_ = cvs;
  }

  const std::vector<TargetedExperiment::CV>& TargetedExperiment::getCVs() const
  {
    return cvs_;
  }

  void TargetedExperiment::addCV(const CV& cv)
  {
    cvs_.push_back(cv);
  }

  void TargetedExperiment::setContacts(const std::vector<Contact>& contacts)
  {
    contacts_ = contacts;
  }

  const std::vector<TargetedExperiment::Contact>& TargetedExperiment::getContacts() const
  {
    return contacts_;
  }

  void TargetedExperiment::addContact(const Contact& contact)
  {
    contacts_.push_back(contact);
  }

  void TargetedExperiment::setPublications(const std::vector<Publication>& publications)
  {
    publications_ = publications;
  }

  const std::vector<TargetedExperiment::Publication>& TargetedExperiment::getPublications() const
  {
    return publications_;
  }

  void TargetedExperiment::addPublication(const Publication& publication)
  {
    publications_.push_back(publication);
  }

  void TargetedExperiment::setTargetCVTerms(const CVTermList& cv_terms)
  {
    targets_ = cv_terms;
  }

  const CVTermList& TargetedExperiment::getTargetCVTerms() const
  {
    return targets_;
  }

  void TargetedExperiment::addTargetCVTerm(const CVTerm& cv_term)
  {
    targets_.addCVTerm(cv_term);
  }

  void TargetedExperiment::setTargetMetaValue(const String &name, const DataValue &value)
  {
    targets_.setMetaValue(name, value);
  }

  void TargetedExperiment::setInstruments(const std::vector<Instrument>& instruments)
  {
    instruments_ = instruments;
  }

  const std::vector<TargetedExperiment::Instrument>& TargetedExperiment::getInstruments() const
  {
    return instruments_;
  }

  void TargetedExperiment::addInstrument(const Instrument& instrument)
  {
    instruments_.push_back(instrument);
  }

  void TargetedExperiment::setSoftware(const std::vector<Software>& software)
  {
    software_ = software;
  }

  const std::vector<Software>& TargetedExperiment::getSoftware() const
  {
    return software_;
  }

  void TargetedExperiment::addSoftware(const Software& software)
  {
    software_.push_back(software);
  }

  void TargetedExperiment::setProteins(const std::vector<Protein>& proteins)
  {
    protein_reference_map_dirty_ = true;
    proteins_ = proteins;
  }

  const std::vector<TargetedExperiment::Protein>& TargetedExperiment::getProteins() const
  {
    return proteins_;
  }

  const TargetedExperiment::Protein& TargetedExperiment::getProteinByRef(const String& ref)
  {
    if (protein_reference_map_dirty_)
    {
      createProteinReferenceMap();
    }
    return *(protein_reference_map_[ref]);
  }

  void TargetedExperiment::addProtein(const Protein& protein)
  {
    protein_reference_map_dirty_ = true;
    proteins_.push_back(protein);
  }

  void TargetedExperiment::setCompounds(const std::vector<Compound>& compounds)
  {
    compounds_ = compounds;
  }

  const std::vector<TargetedExperiment::Compound>& TargetedExperiment::getCompounds() const
  {
    return compounds_;
  }

  void TargetedExperiment::addCompound(const Compound& rhs)
  {
    compounds_.push_back(rhs);
  }

  void TargetedExperiment::setPeptides(const std::vector<Peptide>& peptides)
  {
    peptide_reference_map_dirty_ = true;
    peptides_ = peptides;
  }

  const std::vector<TargetedExperiment::Peptide>& TargetedExperiment::getPeptides() const
  {
    return peptides_;
  }

  const TargetedExperiment::Peptide& TargetedExperiment::getPeptideByRef(const String& ref)
  {
    if (peptide_reference_map_dirty_)
    {
      createPeptideReferenceMap();
    }
    return *(peptide_reference_map_[ref]);
  }

  void TargetedExperiment::addPeptide(const Peptide& rhs)
  {
    peptide_reference_map_dirty_ = true;
    peptides_.push_back(rhs);
  }

  void TargetedExperiment::setTransitions(const std::vector<ReactionMonitoringTransition>& transitions)
  {
    transitions_ = transitions;
  }

  const std::vector<ReactionMonitoringTransition>& TargetedExperiment::getTransitions() const
  {
    return transitions_;
  }

  void TargetedExperiment::addTransition(const ReactionMonitoringTransition& transition)
  {
    transitions_.push_back(transition);
  }

  void TargetedExperiment::setIncludeTargets(const std::vector<IncludeExcludeTarget>& targets)
  {
    include_targets_ = targets;
  }

  const std::vector<IncludeExcludeTarget>& TargetedExperiment::getIncludeTargets() const
  {
    return include_targets_;
  }

  void TargetedExperiment::addIncludeTarget(const IncludeExcludeTarget& target)
  {
    include_targets_.push_back(target);
  }

  void TargetedExperiment::setExcludeTargets(const std::vector<IncludeExcludeTarget>& targets)
  {
    exclude_targets_ = targets;
  }

  const std::vector<IncludeExcludeTarget>& TargetedExperiment::getExcludeTargets() const
  {
    return exclude_targets_;
  }

  void TargetedExperiment::addExcludeTarget(const IncludeExcludeTarget& target)
  {
    exclude_targets_.push_back(target);
  }

  void TargetedExperiment::setSourceFiles(const std::vector<SourceFile>& source_files)
  {
    source_files_ = source_files;
  }

  const std::vector<SourceFile>& TargetedExperiment::getSourceFiles() const
  {
    return source_files_;
  }

  void TargetedExperiment::addSourceFile(const SourceFile& source_file)
  {
    source_files_.push_back(source_file);
  }

  void TargetedExperiment::sortTransitionsByProductMZ()
  {
    std::sort(transitions_.begin(), transitions_.end(), ReactionMonitoringTransition::ProductMZLess() );
  }

  void TargetedExperiment::createProteinReferenceMap() 
  {
    for (Size i = 0; i < getProteins().size(); i++)
    {
      protein_reference_map_[getProteins()[i].id] = &getProteins()[i];
    }

  }

  void TargetedExperiment::createPeptideReferenceMap()
  {
    for (Size i = 0; i < getPeptides().size(); i++)
    {
      peptide_reference_map_[getPeptides()[i].id] = &getPeptides()[i];
    }
  }

}


