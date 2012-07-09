// -*- mode: C++; tab-width: 2; -*-
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
// $Authors: Andreas Bertsch, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/TraMLHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <iostream>

namespace OpenMS
{
  namespace Internal
  {

  TraMLHandler::TraMLHandler(const TargetedExperiment& exp, const String& filename, const String& version, const ProgressLogger& logger)
    : XMLHandler(filename, version),
      logger_(logger),
      exp_(0),
      cexp_(&exp)
  {
    cv_.loadFromOBO("PI",File::find("/CV/psi-ms.obo"));
  }

  TraMLHandler::TraMLHandler(TargetedExperiment& exp, const String& filename, const String& version, const ProgressLogger& logger)
    : XMLHandler(filename, version),
      logger_(logger),
      exp_(&exp),
      cexp_(0)
  {
    cv_.loadFromOBO("PI",File::find("/CV/psi-ms.obo"));
  }  

  TraMLHandler::~TraMLHandler()
  {
  }

  void TraMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
  {

    static const XMLCh* s_type = xercesc::XMLString::transcode("type");
    static const XMLCh* s_value = xercesc::XMLString::transcode("value");
    static const XMLCh* s_name = xercesc::XMLString::transcode("name");

    tag_ = sm_.convert(qname);
    open_tags_.push_back(tag_);

    static std::set<String> tags_to_ignore;
    if ( tags_to_ignore.empty() )
    {
      tags_to_ignore.insert("TraML"); // base node
      tags_to_ignore.insert("ContactList"); // contains only contact sections
      tags_to_ignore.insert("CompoundList"); // contains only compounds
      tags_to_ignore.insert("TransitionList"); // contains only transitions
      tags_to_ignore.insert("ConfigurationList"); // contains only configurations
      tags_to_ignore.insert("cvList"); // contains only CVs
      tags_to_ignore.insert("InstrumentList"); // contains only instruments
      tags_to_ignore.insert("SoftwareList"); // contains only software 
      tags_to_ignore.insert("PublicationList"); // contains only publications
      tags_to_ignore.insert("ProteinList"); // contains only proteins
      tags_to_ignore.insert("SourceFileList"); // contains only source files
      tags_to_ignore.insert("InterpretationList"); // contains only interpretations
      tags_to_ignore.insert("Evidence"); // only cv terms
      tags_to_ignore.insert("ValidationStatus"); // only cv terms
      tags_to_ignore.insert("Sequence"); // only sequence as characters
      tags_to_ignore.insert("Precursor"); // contains only cv terms
      tags_to_ignore.insert("Product"); // contains no attributes
      tags_to_ignore.insert("IntermediateProduct"); // contains no attributes
      tags_to_ignore.insert("TargetIncludeList");
      tags_to_ignore.insert("TargetExcludeList");
      tags_to_ignore.insert("TargetList");
      tags_to_ignore.insert("RetentionTimeList");
    }

    // skip tags where nothing is to do
    if (tags_to_ignore.find(tag_) != tags_to_ignore.end())
    {
      return;
    }

    //determine parent tag
    String parent_tag;
    if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
    String parent_parent_tag;
    if (open_tags_.size()>2) parent_parent_tag = *(open_tags_.end()-3);

    if (tag_ == "cvParam")
    {
      String value, cv_ref, unit_accession, unit_name, unit_cv_ref;
      optionalAttributeAsString_(value, attributes, "value");
      optionalAttributeAsString_(unit_accession, attributes, "unitAccession");
      optionalAttributeAsString_(unit_name, attributes, "unitName");
      optionalAttributeAsString_(unit_cv_ref, attributes, "unitCvRef");
      optionalAttributeAsString_(cv_ref, attributes, "cvRef");
      CVTerm::Unit unit(unit_accession, unit_name, unit_cv_ref);
      CVTerm cv_term(attributeAsString_(attributes, "accession"), attributeAsString_(attributes, "name"), cv_ref, value, unit);

      handleCVParam_(parent_parent_tag, parent_tag, cv_term);
      return;
    }
    else if (tag_=="userParam")
    {
      String type = "";
      optionalAttributeAsString_(type, attributes, s_type);
      String value = "";
      optionalAttributeAsString_(value, attributes, s_value);
      handleUserParam_(parent_parent_tag, parent_tag, attributeAsString_(attributes, s_name), type, value);
    }
    else if (tag_ == "cv")
    {
      exp_->addCV(TargetedExperiment::CV(attributeAsString_(attributes, "id"), attributeAsString_(attributes, "fullName"), attributeAsString_(attributes, "version"), attributeAsString_(attributes, "URI")));
    }
    else if (tag_ == "Contact")
    {
      actual_contact_.id = attributeAsString_(attributes, "id");
    }
    else if (tag_ == "Publication")
    {
      actual_publication_.id = attributeAsString_(attributes, "id");
    }
    else if (tag_ == "Instrument")
    {
      actual_instrument_.id = attributeAsString_(attributes, "id");
    }
    else if (tag_ == "Software")
    {
      actual_software_.setName(attributeAsString_(attributes, "id"));
      actual_software_.setVersion(attributeAsString_(attributes, "version"));
    }
    else if (tag_ == "Protein")
    {
      actual_protein_ = TargetedExperiment::Protein();
      actual_protein_.id = attributeAsString_(attributes, "id");
    }
    else if (tag_ == "Peptide")
    {
      actual_peptide_ = TargetedExperiment::Peptide();
      actual_peptide_.id = attributeAsString_(attributes, "id");
      actual_peptide_.sequence = attributeAsString_(attributes, "sequence");
    }
    else if (tag_ == "Modification")
    {
      TargetedExperiment::Peptide::Modification mod;
      DoubleReal avg_mass_delta(0), mono_mass_delta(0); // zero means no value
      optionalAttributeAsDouble_(avg_mass_delta, attributes, "averageMassDelta");
      optionalAttributeAsDouble_(mono_mass_delta, attributes, "monoisotopicMassDelta");
      mod.avg_mass_delta = avg_mass_delta;
      mod.mono_mass_delta = mono_mass_delta;

      mod.location = attributeAsInt_(attributes, "location");
      actual_peptide_.mods.push_back(mod);
    }
    else if (tag_ == "Compound")
    {
      actual_compound_ = TargetedExperiment::Compound();
      actual_compound_.id = attributeAsString_(attributes, "id");
    }
    else if (tag_ == "Prediction")
    {
      actual_prediction_.software_ref = attributeAsString_(attributes, "softwareRef");
      String contact_ref;
      if (optionalAttributeAsString_(contact_ref, attributes, "contactRef"))
      {
        actual_prediction_.contact_ref = contact_ref;
      }
    }
    else if (tag_ == "RetentionTime")
    {
      actual_rt_ = TargetedExperiment::RetentionTime();
      String software_ref;
      if (optionalAttributeAsString_(software_ref, attributes, "softwareRef"))
      {
        actual_rt_.software_ref = software_ref;
      }
    }
    else if (tag_ == "Transition")
    {
      actual_transition_ = ReactionMonitoringTransition();
      String id;
      if (optionalAttributeAsString_(id, attributes, "id"))
      {
        actual_transition_.setName(id);
      }
      String peptide_ref;
      if (optionalAttributeAsString_(peptide_ref, attributes, "peptideRef"))
      {
        actual_transition_.setPeptideRef(peptide_ref);
      }
      String compound_ref;
      if (optionalAttributeAsString_(compound_ref, attributes, "compoundRef"))
      {
        actual_transition_.setCompoundRef(compound_ref);
      }
    }
    else if (tag_ == "Interpretation")
    {
      String primary;
      if (optionalAttributeAsString_(primary, attributes, "primary"))
      {
        actual_interpretation_.setMetaValue("primary", primary);
      }
    }
    else if (tag_ == "Configuration")
    {
      actual_configuration_.instrument_ref = attributeAsString_(attributes, "instrumentRef");
      String contact_ref;
      if (optionalAttributeAsString_(contact_ref, attributes, "contactRef"))
      {
        actual_configuration_.contact_ref = contact_ref;
      }
    }
    else if (tag_ == "SourceFile")
    {
      actual_sourcefile_.setNativeIDType(attributeAsString_(attributes, "id"));
      actual_sourcefile_.setNameOfFile(attributeAsString_(attributes, "name"));
      actual_sourcefile_.setPathToFile(attributeAsString_(attributes, "location"));
    }
    else if (tag_ == "ProteinRef")
    {
      actual_peptide_.protein_refs.push_back(attributeAsString_(attributes, "ref"));
    }
    else if (tag_ == "Target")
    {
      actual_target_ = IncludeExcludeTarget();
      String id;
      if (optionalAttributeAsString_(id, attributes, "id"))
      {
        actual_target_.setName(id);
      }
      String peptide_ref;
      if (optionalAttributeAsString_(peptide_ref, attributes, "peptideRef"))
      {
        actual_target_.setPeptideRef(peptide_ref);
      }
      String compound_ref;
      if (optionalAttributeAsString_(compound_ref, attributes, "compoundRef"))
      {
        actual_target_.setCompoundRef(compound_ref);
      }
    }
    else 
    {
      error(LOAD, "TraMLHandler: unknown tag opening: '" + tag_ + "'");
    }
    return;
  }

  void TraMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
  {
    if (open_tags_.back() == "Sequence")
    {
      String protein_sequence = sm_.convert(chars);
      actual_protein_.sequence = protein_sequence;
      return;
    }
    return;
  }

  void TraMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    tag_ = sm_.convert(qname);

    //determine parent tag
    String parent_tag;
    if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
    String parent_parent_tag;
    if (open_tags_.size()>2) parent_parent_tag = *(open_tags_.end()-3);

    open_tags_.pop_back();

    static std::set<String> tags_to_ignore;
    if ( tags_to_ignore.empty() )
    {
      tags_to_ignore.insert("TraML"); // base node
      tags_to_ignore.insert("ContactList"); // contains only contact sections
      tags_to_ignore.insert("CompoundList"); // contains only compounds
      tags_to_ignore.insert("TransitionList"); // contains only transitions
      tags_to_ignore.insert("ConfigurationList"); // contains only configurations
      tags_to_ignore.insert("cvList"); // contains only CVs
      tags_to_ignore.insert("InstrumentList"); // contains only instruments
      tags_to_ignore.insert("SoftwareList"); // contains only software 
      tags_to_ignore.insert("PublicationList"); // contains only publications
      tags_to_ignore.insert("ProteinList"); // contains only proteins
      tags_to_ignore.insert("SourceFileList"); // contains only source files
      tags_to_ignore.insert("InterpretationList"); // contains only interpretations
      tags_to_ignore.insert("Evidence"); // only cv terms
      tags_to_ignore.insert("cvParam"); // already handled
      tags_to_ignore.insert("userParam"); // already handled
      tags_to_ignore.insert("cv"); // already handled
      tags_to_ignore.insert("Sequence"); // already handled in characters
      tags_to_ignore.insert("Precursor"); // contains only cv terms
      tags_to_ignore.insert("RetentionTimeList");
      tags_to_ignore.insert("TargetList");
      tags_to_ignore.insert("TargetIncludeList");
      tags_to_ignore.insert("TargetExcludeList");  
      tags_to_ignore.insert("ProteinRef");
      tags_to_ignore.insert("Modification");
      tags_to_ignore.insert("TargetList");
    }

    // skip tags where nothing is to do
    if (tags_to_ignore.find(tag_) != tags_to_ignore.end())
    {
      return;
    }
    else if (tag_ == "Contact")
    {
      exp_->addContact(actual_contact_);
      actual_contact_ = TargetedExperiment::Contact();
    }
    else if (tag_ == "Instrument")
    {
      exp_->addInstrument(actual_instrument_);
      actual_instrument_ = TargetedExperiment::Instrument();
    }
    else if (tag_ == "Publication")
    {
      exp_->addPublication(actual_publication_);
      actual_publication_ = TargetedExperiment::Publication();
    }
    else if (tag_ == "Software")
    {
      exp_->addSoftware(actual_software_);
      actual_software_ = Software();
    }
    else if (tag_ == "Protein")
    {
      exp_->addProtein(actual_protein_);
    }
    else if (tag_ == "RetentionTime")
    {
      if (parent_parent_tag == "Peptide")
      {
        actual_peptide_.rts.push_back(actual_rt_);
        actual_rt_ = TargetedExperiment::RetentionTime();
      }
      else if (parent_parent_tag == "Compound")
      {
        actual_compound_.rts.push_back(actual_rt_);
        actual_rt_ = TargetedExperiment::RetentionTime();
      }
      else if (parent_tag == "Target")
      {
        actual_target_.setRetentionTime(actual_rt_);
        actual_rt_ = TargetedExperiment::RetentionTime();
      }
      else if (parent_tag == "Transition")
      {
        actual_transition_.setRetentionTime(actual_rt_);
        actual_rt_ = TargetedExperiment::RetentionTime();
      }
      else 
      {
        error(LOAD, "TraMLHandler: tag 'RetentionTime' not allowed at parent tag '" + parent_tag + "', ignoring!");
      }
    }
    else if (tag_ == "Peptide")
    {
      exp_->addPeptide(actual_peptide_);
      actual_peptide_ = TargetedExperiment::Peptide();
    }
    else if (tag_ == "Compound")
    {
      exp_->addCompound(actual_compound_);
      actual_compound_ = TargetedExperiment::Compound();
    }
    else if (tag_ == "Transition")
    {
      exp_->addTransition(actual_transition_);
      actual_transition_ = ReactionMonitoringTransition();
    }
    else if (tag_ == "Product")
    {
      actual_transition_.setProduct(actual_product_);
      actual_product_ = ReactionMonitoringTransition::Product();
    }
    else if (tag_ == "IntermediateProduct")
    {
      actual_transition_.addIntermediateProduct(actual_product_);
      actual_product_ = ReactionMonitoringTransition::Product();
    }
    else if (tag_ == "Interpretation")
    {
      actual_product_.addInterpretation(actual_interpretation_);
      actual_interpretation_ = CVTermList();
    }
    else if (tag_ == "Prediction")
    {
      actual_transition_.setPrediction(actual_prediction_);
      actual_prediction_ = TargetedExperiment::Prediction();
    }
    else  if (tag_ == "Configuration")
    {
      if (parent_parent_tag == "IntermediateProduct" || parent_parent_tag == "Product")
      {
        actual_product_.addConfiguration(actual_configuration_);
        actual_configuration_ = TargetedExperimentHelper::Configuration();
      }
      else if (parent_parent_tag == "Target")
      {
        actual_target_.addConfiguration(actual_configuration_);
        actual_configuration_ = TargetedExperimentHelper::Configuration();
      }
      else 
      {
        error(LOAD, "TraMLHandler: tag 'Configuration' not allowed at parent tag '" + parent_tag + "', ignoring!");
      }
    }
    else if (tag_ == "ValidationStatus")
    {
      actual_configuration_.validations.push_back(actual_validation_);
      actual_validation_ = CVTermList();
    }
    else if (tag_ == "SourceFile")
    {
      exp_->addSourceFile(actual_sourcefile_);
      actual_sourcefile_ = SourceFile();
    }
    else if (tag_ == "Target")
    {
      if (parent_tag == "TargetIncludeList")
      {
        exp_->addIncludeTarget(actual_target_);
        actual_target_ = IncludeExcludeTarget();
      }
      else if (parent_tag == "TargetExcludeList")
      {
        exp_->addExcludeTarget(actual_target_);
        actual_target_ = IncludeExcludeTarget();
      }
      else
      {
        error(LOAD, "TraMLHandler: tag 'Target' not allowed at parent tag '" + parent_tag + "', ignoring!");
      }
    }
    else
    {
      error(LOAD, "TraMLHandler: unknown tag closing: '" + tag_ + "'");
    }
    return;
  }
  
  void TraMLHandler::writeTo(std::ostream& os)
  {
    const TargetedExperiment& exp = *(cexp_);
    logger_.startProgress(0,exp.getTransitions().size(),"storing TraML file");
    int progress = 0;

    os  << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n";
    os  << "<TraML version=\"1.0.0\" xmlns=\"http://psi.hupo.org/ms/traml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/traml TraML1.0.0.xsd\">" << "\n";

    //--------------------------------------------------------------------------------------------
    // CV list
    //--------------------------------------------------------------------------------------------
    os  << "  <cvList>" << "\n";

    if (exp.getCVs().size() == 0)
    {
      os  << "    <cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"unknown\" URI=\"http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\"/>" << "\n"
          << "    <cv id=\"UO\" fullName=\"Unit Ontology\" version=\"unknown\" URI=\"http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo\"/>" << "\n";
    }
    else
    {
      for (std::vector<TargetedExperiment::CV>::const_iterator it = exp.getCVs().begin(); it != exp.getCVs().end(); ++it)
      {
        os << "    <cv id=\"" << it->id << "\" fullName=\"" << it->fullname << "\" version=\"" << it->version << "\" URI=\"" << it->URI << "\"/>" << "\n";
      }
    }
    os  << "  </cvList>" << "\n";

    // source file list
    if (exp.getSourceFiles().size() > 0)
    {
      os << "  <SourceFileList>" << "\n";
      for (std::vector<SourceFile>::const_iterator it = exp.getSourceFiles().begin(); it != exp.getSourceFiles().end(); ++it)
      {
        os << "    <SourceFile id=\"" 
          << it->getNativeIDType() << "\" name=\"" 
          << it->getNameOfFile() << "\" location=\"" 
          << it->getPathToFile() << "\">" 
          << "\n";
        writeCVParams_(os, *it, 3);
        writeUserParam_(os, (MetaInfoInterface)*it, 3);
        os << "    </SourceFile>" << "\n";
      }
      os << "  </SourceFileList>" << "\n";
    }

    // contact list
    if (exp.getContacts().size() > 0)
    {
      os << "  <ContactList>" << "\n";
      for (std::vector<TargetedExperiment::Contact>::const_iterator it = exp.getContacts().begin(); it != exp.getContacts().end(); ++it)
      {
        os << "    <Contact id=\"" << it->id << "\">" << "\n";
        writeCVParams_(os, *it, 3);
        writeUserParam_(os, (MetaInfoInterface)*it, 3);
        os << "    </Contact>" << "\n";
      }
      os << "  </ContactList>" << "\n";
    }

    // publication list
    if (exp.getPublications().size() > 0)
    {
      os << "  <PublicationList>"  << "\n";
      for (std::vector<TargetedExperiment::Publication>::const_iterator it = exp.getPublications().begin(); it != exp.getPublications().end(); ++it)
      {  
        os << "    <Publication id=\"" << it->id << "\">" << "\n";
        writeCVParams_(os, *it, 3);
        writeUserParam_(os, (MetaInfoInterface)*it, 3);
        os << "    </Publication>" << "\n";
      }
      os << "  </PublicationList>" << "\n";
    }

    // instrument list
    if (exp.getInstruments().size() > 0)
    {
      os << "  <InstrumentList>" << "\n";
      for (std::vector<TargetedExperiment::Instrument>::const_iterator it = exp.getInstruments().begin(); it != exp.getInstruments().end(); ++it)
      {
        os << "    <Instrument id=\"" << it->id << "\">" << "\n";
        writeCVParams_(os, *it, 3);
        writeUserParam_(os, (MetaInfoInterface)*it, 3);
        os << "    </Instrument>" << "\n";
      }
      os << "  </InstrumentList>" << "\n";
    }

    // software list
    if (exp.getSoftware().size() > 0)
    {
      os << "  <SoftwareList>" << "\n";
      for (std::vector<Software>::const_iterator it = exp.getSoftware().begin(); it != exp.getSoftware().end(); ++it)
      {
        os << "    <Software id=\"" << it->getName() << "\" version=\"" << it->getVersion() << "\">" << "\n";
        writeCVParams_(os, (CVTermList)*it, 3);
        writeUserParam_(os, (MetaInfoInterface)*it, 3);
        os << "    </Software>" << "\n";
      }
      os << "  </SoftwareList>" << "\n";
    }

    //--------------------------------------------------------------------------------------------
    // protein list
    //--------------------------------------------------------------------------------------------
    if (exp.getProteins().size() > 0)
    {
      os << "  <ProteinList>" << "\n";
      for (std::vector<TargetedExperiment::Protein>::const_iterator it = exp.getProteins().begin(); it != exp.getProteins().end(); ++it)
      {
        os << "    <Protein id=\"" << it->id << "\">" << "\n";
        writeCVParams_(os, (CVTermList)*it, 3);
        writeUserParam_(os, (MetaInfoInterface)*it, 3);
        os << "      <Sequence>" << it->sequence << "</Sequence>" << "\n";
        os << "    </Protein>" << "\n";
      }
      os << "  </ProteinList>" << "\n";
    }

    //--------------------------------------------------------------------------------------------
    // compound list
    //--------------------------------------------------------------------------------------------
    if (exp.getCompounds().size()  + exp.getPeptides().size() > 0)
    {
      os << "  <CompoundList>" << "\n";
      std::vector<TargetedExperiment::Peptide> exp_peptides = exp.getPeptides();
      
      for (std::vector<TargetedExperiment::Peptide>::const_iterator it = exp_peptides.begin(); it != exp_peptides.end(); ++it)
      {
        os << "    <Peptide id=\"" << it->id << "\" sequence=\"" << it->sequence << "\">" << "\n";
        if ( it->getChargeState() != -1 )
        {
          os << "      <cvParam cvRef=\"MS\" accession=\"MS:1000041\" name=\"charge state\" value=\"" <<  it->getChargeState() << "\"/>\n";
        }
        if ( it->getPeptideGroupLabel() != "" )
        {
          os << "      <cvParam cvRef=\"MS\" accession=\"MS:1000893\" name=\"peptide group label\" value=\"" <<  it->getPeptideGroupLabel() << "\"/>\n";
        }
        writeCVParams_(os, (CVTermList)*it, 3);
        writeUserParam_(os, (MetaInfoInterface)*it, 3);

        for (std::vector<String>::const_iterator rit = it->protein_refs.begin(); rit != it->protein_refs.end(); ++rit)
        {
          os << "      <ProteinRef ref=\"" << *rit << "\"/>" << "\n";
        }

        if (it->mods.size() > 0)
        {  
          for (std::vector<TargetedExperiment::Peptide::Modification>::const_iterator mit = it->mods.begin(); mit != it->mods.end(); ++mit)
          {
            os << "      <Modification";  
            os << " location=\"" << mit->location << "\""; // location is required
            if (mit->mono_mass_delta != 0)
            {
              os << " monoisotopicMassDelta=\"" << mit->mono_mass_delta << "\"";
            }
            if (mit->avg_mass_delta != 0)
            {
              os << " averageMassDelta=\"" << mit->avg_mass_delta << "\"";
            }
            os << ">\n";  
            writeCVParams_(os, (CVTermList)*mit, 4);
            writeUserParam_(os, (MetaInfoInterface)*mit, 4);
            os << "      </Modification>\n";  
          }
        }
      
        if (it->rts.size() > 0)
        {  
          os << "      <RetentionTimeList>\n";  
          for (std::vector<TargetedExperiment::RetentionTime>::const_iterator rit = it->rts.begin(); rit != it->rts.end(); ++rit)
          {
            os << "        <RetentionTime";
            if (rit->software_ref != "")
            {
              os << " softwareRef=\"" << rit->software_ref << "\"";
            }
            os << ">" << "\n";
            writeCVParams_(os, (CVTermList)*rit, 5);
            writeUserParam_(os, (MetaInfoInterface)*rit, 5);
            os << "        </RetentionTime>" << "\n";
          }
          os << "      </RetentionTimeList>\n";  
        }
  
        if (!it->evidence.empty())
        {
          os << "      <Evidence>" << "\n";
          writeCVParams_(os, it->evidence, 4);
          writeUserParam_(os, (MetaInfoInterface)it->evidence, 4);
          os << "      </Evidence>" << "\n";
        }
        os << "    </Peptide>" << "\n";
      }

      for (std::vector<TargetedExperiment::Compound>::const_iterator it = exp.getCompounds().begin(); it != exp.getCompounds().end(); ++it)
      {
        os << "    <Compound id=\"" << it->id << "\">" << "\n";
        writeCVParams_(os, (CVTermList)*it, 3);
        writeUserParam_(os, (MetaInfoInterface)*it, 3);

        if (it->rts.size() > 0)
        {  
          os << "      <RetentionTimeList>\n";  
          for (std::vector<TargetedExperiment::RetentionTime>::const_iterator rit = it->rts.begin(); rit != it->rts.end(); ++rit)
            {
              os << "        <RetentionTime";
              if (rit->software_ref != "")
              {
               os << " softwareRef=\"" << rit->software_ref << "\"";
              }
              os << ">" << "\n";
              writeCVParams_(os, (CVTermList)*rit, 5);
            writeUserParam_(os, (MetaInfoInterface)*rit, 5);
              os << "        </RetentionTime>" << "\n";
            }
          os << "      </RetentionTimeList>\n";  
        }
        os << "    </Compound>" << "\n";
      }

      os << "  </CompoundList>" << "\n";
    }

    //--------------------------------------------------------------------------------------------
    // transition list
    //--------------------------------------------------------------------------------------------
    if (exp.getTransitions().size() > 0)
    {
      os << "  <TransitionList>" << "\n";
      for (std::vector<ReactionMonitoringTransition>::const_iterator it = exp.getTransitions().begin(); it != exp.getTransitions().end(); ++it)
      {
        logger_.setProgress(progress++);
        os << "    <Transition";
        os << " id=\"" << it->getName() << "\"";

        if (it->getPeptideRef() != "")
        {
          os << " peptideRef=\"" << it->getPeptideRef() << "\"";
        }

        if (it->getCompoundRef() != "")
        {
          os << " compoundRef=\"" << it->getCompoundRef() << "\"";
        }
        os << ">" << "\n";


        if ( it->getLibraryIntensity() > -100 )
        {
          os << "      <cvParam cvRef=\"MS\" accession=\"MS:1001226\" name=\"product ion intensity\" value=\"" <<  it->getLibraryIntensity() << "\"/>\n";
        }
        if ( it->getDecoyTransitionType() != ReactionMonitoringTransition::UNKNOWN )
        {
          if ( it->getDecoyTransitionType() == ReactionMonitoringTransition::TARGET )
          {
            os << "      <cvParam cvRef=\"MS\" accession=\"MS:1002007\" name=\"target SRM transition\"/>\n";
          }
          else if ( it->getDecoyTransitionType() == ReactionMonitoringTransition::DECOY )
          {
            os << "      <cvParam cvRef=\"MS\" accession=\"MS:1002008\" name=\"decoy SRM transition\"/>\n";
          }
        }


        writeCVParams_(os, (CVTermList)*it, 3);
        writeUserParam_(os, (MetaInfoInterface)*it, 3);

        // Precursor is required
        os << "      <Precursor>" << "\n"; 
        os << "        <cvParam cvRef=\"MS\" accession=\"MS:1000827\" name=\"isolation window target m/z\" value=\"" << precisionWrapper(it->getPrecursorMZ()) << "\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n";
        writeCVParams_(os, it->getPrecursorCVTermList(), 4);
        writeUserParam_(os, (MetaInfoInterface)it->getPrecursorCVTermList(), 4);
        os << "      </Precursor>" << "\n";
      
        for (ProductListType::const_iterator prod_it = it->getIntermediateProducts().begin(); prod_it != it->getIntermediateProducts().end(); prod_it++)
        {
          os << "      <IntermediateProduct>" << "\n";
          write_product_ (os, prod_it);
          os << "      </IntermediateProduct>" << "\n";
        }

        // Product is required
        os << "      <Product>" << "\n";
        ProductListType dummy_vect;
        dummy_vect.push_back(it->getProduct());
        write_product_ (os, dummy_vect.begin() );
        os << "      </Product>" << "\n";

        const IncludeExcludeTarget::RetentionTime* rit = &it->getRetentionTime();
        if (!rit->getCVTerms().empty())
        {
          os << "      <RetentionTime";
          if (rit->software_ref != "")
          {
            os << " softwareRef=\"" << rit->software_ref << "\"";
          }
          os << ">" << "\n";
          writeCVParams_(os, (CVTermList)*rit, 4);
          writeUserParam_(os, (MetaInfoInterface)*rit, 4);
          os << "      </RetentionTime>" << "\n";
        }

        if (!it->getPrediction().empty())
        {
          os << "      <Prediction softwareRef=\"" << it->getPrediction().software_ref << "\"";
          if (!it->getPrediction().contact_ref.empty())
          {
            os << " contactRef=\"" << it->getPrediction().contact_ref << "\"";
          }
          os << ">" << "\n";
          writeCVParams_(os, it->getPrediction(), 4);
          writeUserParam_(os, (MetaInfoInterface)it->getPrediction(), 4);
          os << "      </Prediction>" << "\n";
        }
        
        os << "    </Transition>" << "\n";
      }
      os << "  </TransitionList>" << "\n";
    
    }

    if (!exp.getIncludeTargets().empty() || !exp.getExcludeTargets().empty() )
    {
      os << "  <TargetList>" << "\n";
      writeCVParams_(os, exp.getTargetCVTerms(), 2);
      writeUserParam_(os, (MetaInfoInterface)exp.getTargetCVTerms(), 2);

      if (!exp.getIncludeTargets().empty())
      {
        os << "    <TargetIncludeList>" << "\n";
        for (std::vector<IncludeExcludeTarget>::const_iterator it = exp.getIncludeTargets().begin(); it != exp.getIncludeTargets().end(); ++it)
        {
          write_target_(os, it);
        }
        os << "    </TargetIncludeList>" << "\n";
      }

      if (!exp.getExcludeTargets().empty() )
      {
        os << "    <TargetExcludeList>" << "\n";
        for (std::vector<IncludeExcludeTarget>::const_iterator it = exp.getExcludeTargets().begin(); it != exp.getExcludeTargets().end(); ++it)
        {
          write_target_(os, it);
        }
        os << "    </TargetExcludeList>" << "\n";
      }

      os << "  </TargetList>" << "\n";
    }

    os << "</TraML>" << "\n";
    logger_.endProgress();
    return;
  }

  void TraMLHandler::write_target_ (std::ostream& os, const std::vector<IncludeExcludeTarget>::const_iterator& it) const
  {
      os << "      <Target id=\"" << it->getName() << "\"";
      if (!it->getPeptideRef().empty())
      {
        os << " peptideRef=\"" << it->getPeptideRef() << "\"";
      }
      if (!it->getCompoundRef().empty())
      {
        os << " compoundRef=\"" << it->getCompoundRef() << "\"";
      }
      os << ">\n";
      os << "        <Precursor>\n";
      writeCVParams_(os, it->getPrecursorCVTermList(), 5);
      writeUserParam_(os, (MetaInfoInterface)it->getPrecursorCVTermList(), 5);
      os << "        </Precursor>\n";

      const IncludeExcludeTarget::RetentionTime* rit = &it->getRetentionTime();
      if (!rit->getCVTerms().empty())
      {
        os << "        <RetentionTime";
        if (rit->software_ref != "")
        {
          os << " softwareRef=\"" << rit->software_ref << "\"";
        }
        os << ">" << "\n";
        writeCVParams_(os, (CVTermList)*rit, 5);
        writeUserParam_(os, (MetaInfoInterface)*rit, 5);
        os << "        </RetentionTime>" << "\n";
      }

      if (!it->getConfigurations().empty())
      {
        os << "        <ConfigurationList>\n";
        for (std::vector<TargetedExperimentHelper::Configuration>::const_iterator config_it = it->getConfigurations().begin(); config_it != it->getConfigurations().end(); ++config_it)
        {
          write_configuration_ (os, config_it);
        }
        os << "        </ConfigurationList>\n";
      }
      os << "      </Target>" << "\n";

  }

  void TraMLHandler::write_product_ (std::ostream& os, const std::vector<ReactionMonitoringTransition::Product>::const_iterator& prod_it) const
  {
    writeCVParams_(os, *prod_it, 4);
    writeUserParam_(os, (MetaInfoInterface)*prod_it, 4);

    if (!prod_it->getInterpretationList().empty())
    {
      os << "        <InterpretationList>" << "\n";
      for (std::vector<CVTermList>::const_iterator inter_it = prod_it->getInterpretationList().begin(); inter_it != prod_it->getInterpretationList().end(); inter_it++)
      {
        os << "          <Interpretation>" << "\n";
        writeCVParams_(os, *inter_it, 6);
        writeUserParam_(os, (MetaInfoInterface)*inter_it, 6);
        os << "          </Interpretation>" << "\n";
      }
      os << "        </InterpretationList>" << "\n";
    }
    if (!prod_it->getConfigurationList().empty())
    {
      os << "        <ConfigurationList>" << "\n";
      for (ConfigurationListType::const_iterator config_it = prod_it->getConfigurationList().begin(); config_it != prod_it->getConfigurationList().end(); config_it++)
      {
        write_configuration_ (os, config_it);
      }
      os << "        </ConfigurationList>" << "\n";
    }
  }

  void TraMLHandler::write_configuration_ (std::ostream& os, const std::vector<ReactionMonitoringTransition::Configuration>::const_iterator& cit) const
  {
    os << "          <Configuration instrumentRef=\"" << cit->instrument_ref << "\"";
    if (cit->contact_ref != "")
    {
      os << " contactRef=\"" << cit->contact_ref << "\"";
    }
    os << ">" << "\n";

    writeCVParams_(os, (CVTermList)*cit, 6);
    writeUserParam_(os, (MetaInfoInterface)*cit, 6);
    if (cit->validations.size() != 0)
    {
      for (std::vector<CVTermList>::const_iterator iit = cit->validations.begin(); iit != cit->validations.end(); ++iit)
      {
        if (!iit->empty())
        {
          os << "            <ValidationStatus>" << "\n";
          writeCVParams_(os, *iit, 7);
          writeUserParam_(os, (MetaInfoInterface)*iit, 7);
          os << "            </ValidationStatus>" << "\n";
        }
      }
    }

    os << "          </Configuration>" << "\n";
  }

  void TraMLHandler::writeCVParams_(std::ostream& os, const CVTermList& cv_terms, UInt indent) const
  {
    for (Map<String, std::vector<CVTerm> >::const_iterator it = cv_terms.getCVTerms().begin(); it != cv_terms.getCVTerms().end(); ++it)
    {
      for (std::vector<CVTerm>::const_iterator cit = it->second.begin(); cit != it->second.end(); ++cit)
      {
        os << String(2 * indent, ' ') << "<cvParam cvRef=\"" << cit->getCVIdentifierRef() << "\" accession=\"" << cit->getAccession() << "\" name=\"" << cit->getName() << "\"";
        if (cit->hasValue() && !cit->getValue().isEmpty() && !cit->getValue().toString().empty())
        {
          os << " value=\"" << cit->getValue().toString() << "\"";
        }

        if (cit->hasUnit())
        {
          os << " unitCvRef=\"" << cit->getUnit().cv_ref << "\" unitAccession=\"" << cit->getUnit().accession << "\" unitName=\"" << cit->getUnit().name << "\"";
        }
        os << "/>" << "\n";
      }
    }
  }

  void TraMLHandler::handleCVParam_(const String& parent_parent_tag, const String& parent_tag, const CVTerm& cv_term)
  {
      //Error checks of CV values
      String accession = cv_term.getAccession();
      if (cv_.exists(accession))
      {
        const ControlledVocabulary::CVTerm& term = cv_.getTerm(accession);
        //obsolete CV terms
        if (term.obsolete)
        {
          warning(LOAD, String("Obsolete CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "'.");
        }
        //check if term name and parsed name match
        String parsed_name = cv_term.getName();
        parsed_name.trim();
        String correct_name = term.name;
        correct_name.trim();
        if (parsed_name!=correct_name)
        {
          warning(LOAD, String("Name of CV term not correct: '") + term.id + " - " + parsed_name + "' should be '" + correct_name + "'");
        }
        if (term.obsolete)
        {
          warning(LOAD, String("Obsolete CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "'.");
        //values used in wrong places and wrong value types
        String value = cv_term.getValue().toString();
        if (value != "")
        {
          if (term.xref_type==ControlledVocabulary::CVTerm::NONE)
          {
            //Quality CV does not state value type :(
            if (!accession.hasPrefix("PATO:"))
            {
              warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must not have a value. The value is '" + value + "'.");
            }
          }
          else
          {
            switch (term.xref_type)
            {
              //string value can be anything
              case ControlledVocabulary::CVTerm::XSD_STRING:
                break;
              //int value => try casting
              case ControlledVocabulary::CVTerm::XSD_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_NEGATIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_POSITIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_NON_NEGATIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_NON_POSITIVE_INTEGER:
                try
                {
                  value.toInt();
                }
                catch(Exception::ConversionError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must have an integer value. The value is '" + value + "'.");
                  return;
                }
                break;
              //double value => try casting
              case ControlledVocabulary::CVTerm::XSD_DECIMAL:
                try
                {
                  value.toDouble();
                }
                catch(Exception::ConversionError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must have a floating-point value. The value is '" + value + "'.");
                  return;
                }
                break;
              //date string => try conversion
              case ControlledVocabulary::CVTerm::XSD_DATE:
                try
                {
                  DateTime tmp;
                  tmp.set(value);
                }
                catch(Exception::ParseError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must be a valid date. The value is '" + value + "'.");
                  return;
                }
                break;
              default:
                warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' has the unknown value type '" + ControlledVocabulary::CVTerm::getXRefTypeName(term.xref_type) + "'.");
                break;
            }
          }
        }
        //no value, although there should be a numerical value
        else if (term.xref_type!=ControlledVocabulary::CVTerm::NONE && term.xref_type!=ControlledVocabulary::CVTerm::XSD_STRING)
        {
          warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' should have a numerical value. The value is '" + value + "'.");
          return;
        }
      }
    }


    // now handle the CVTerm and add it to the object
    if (parent_tag == "Software")
    {
      actual_software_.addCVTerm(cv_term);
    }
    else if (parent_tag == "Publication")
    {
      actual_publication_.addCVTerm(cv_term);
    } 
    else if (parent_tag == "Instrument")
    {
      actual_instrument_.addCVTerm(cv_term);
    } 
    else if (parent_tag == "Contact")
    {
      actual_contact_.addCVTerm(cv_term);
    }
    else if (parent_tag == "RetentionTime")
    {
      actual_rt_.addCVTerm(cv_term);
    }
    else if (parent_tag == "Evidence")
    {
      actual_peptide_.evidence.addCVTerm(cv_term);
    }
    else if (parent_tag == "Peptide")
    {
      if (cv_term.getAccession() == "MS:1000041")
      {
        actual_peptide_.setChargeState(cv_term.getValue().toString().toInt());
      }
      else if (cv_term.getAccession() == "MS:1000893")
      {
        actual_peptide_.setPeptideGroupLabel(cv_term.getValue().toString());
      }
      else
      {
        actual_peptide_.addCVTerm(cv_term);
      }
    }
    else if (parent_tag == "Modification")
    {
      actual_peptide_.mods.back().addCVTerm(cv_term);
    }
    else if (parent_tag == "Compound")
    {
      actual_compound_.addCVTerm(cv_term);
    }
    else if (parent_tag == "Protein")
    {
      actual_protein_.addCVTerm(cv_term);
    }
    else if (parent_tag == "Configuration")
    {
      actual_configuration_.addCVTerm(cv_term);
    }
    else if (parent_tag == "Prediction")
    {
      actual_prediction_.addCVTerm(cv_term);
    }
    else if (parent_tag == "Interpretation")
    {
      actual_interpretation_.addCVTerm(cv_term);
    }
    else if (parent_tag == "ValidationStatus")
    {
      actual_validation_.addCVTerm(cv_term);
    }
    else if (parent_tag == "TargetList")
    {
      exp_->addTargetCVTerm(cv_term);
    }
    else if (parent_tag == "Target")
    {
      actual_target_.addCVTerm(cv_term);
    }
    else if (parent_tag == "Precursor")
    {
      if (parent_parent_tag == "Transition")
      {
        // handle specific CV terms of Transition, currently these are 
        // id: MS:1000827 name: isolation window target m/z
        if (cv_term.getAccession() == "MS:1000827")
        {
          actual_transition_.setPrecursorMZ(cv_term.getValue().toString().toDouble());
        }
        else 
        {
          actual_transition_.addPrecursorCVTerm(cv_term);
        }
      }
      if (parent_parent_tag == "Target")
      {
        actual_target_.addPrecursorCVTerm(cv_term);
      }
    }
    else if (parent_tag == "IntermediateProduct")
    {
      actual_product_.addCVTerm(cv_term);
    }
    else if (parent_tag == "Product")
    {
      actual_product_.addCVTerm(cv_term);
    }
    else if (parent_tag == "SourceFile")
    {
      // TODO handle checksum type...
      actual_sourcefile_.addCVTerm(cv_term);
    } 
    else if (parent_tag == "Transition")
    {
      // handle specific CV terms of Transition, currently these are 
      // id: MS:1002007 name: target SRM transition
      // id: MS:1002008 name: decoy SRM transition
      //
      // id: MS:1000905 (percent of base peak times 100) or MS:1001226 (product ion intensity)
      if (cv_term.getAccession() == "MS:1002007")
      {
        actual_transition_.setDecoyTransitionType(ReactionMonitoringTransition::TARGET);
      }
      else if (cv_term.getAccession() == "MS:1002008")
      {
        actual_transition_.setDecoyTransitionType(ReactionMonitoringTransition::DECOY);
      }
      else if (cv_term.getAccession() == "MS:1001226")
      {
        actual_transition_.setLibraryIntensity(cv_term.getValue().toString().toDouble());
      }
      else if (cv_term.getAccession() == "MS:1000905")
      {
        actual_transition_.setLibraryIntensity(cv_term.getValue().toString().toDouble());
      }
      else
      {
        actual_transition_.addCVTerm(cv_term);
      }
    }
    else 
    {
      warning(LOAD, String("The CV term '" + cv_term.getAccession() + "' - '" + cv_term.getName() + "' used in tag '" + parent_tag + "' could not be handled, ignoring it!"));
    }
    return;

  }

  void TraMLHandler::handleUserParam_(const String& parent_parent_tag, const String& parent_tag, const String& name, const String& type, const String& value)
  {
    //create a DataValue that contains the data in the right type
    DataValue data_value;
    //float type
    if (type=="xsd:double" || type=="xsd:float")
    {
      data_value = DataValue(value.toDouble());
    }
    //integer type
    else if (type=="xsd:byte" || type=="xsd:decimal" || type=="xsd:int" || type=="xsd:integer" || type=="xsd:long" || type=="xsd:negativeInteger" || type=="xsd:nonNegativeInteger" || type=="xsd:nonPositiveInteger" || type=="xsd:positiveInteger" || type=="xsd:short" || type=="xsd:unsignedByte" || type=="xsd:unsignedInt" || type=="xsd:unsignedLong" || type=="xsd:unsignedShort")
    {
      data_value = DataValue(value.toInt());
    }
    //everything else is treated as a string
    else
    {
      data_value = DataValue(value);
    }
    
    //find the right MetaInfoInterface
    if (parent_tag == "Software")
    {
      actual_software_.setMetaValue(name, data_value);
    }
    else if (parent_tag == "Publication")
    {
      actual_publication_.setMetaValue(name, data_value);
    } 
    else if (parent_tag == "Instrument")
    {
      actual_instrument_.setMetaValue(name, data_value);
    } 
    else if (parent_tag == "Contact")
    {
      actual_contact_.setMetaValue(name, data_value);
    }
    else if (parent_tag == "RetentionTime")
    {
      actual_rt_.setMetaValue(name, data_value);
    }
    else if (parent_tag == "Evidence")
    {
      actual_peptide_.evidence.setMetaValue(name, data_value);
    }
    else if (parent_tag == "Peptide")
    {
      actual_peptide_.setMetaValue(name, data_value);
    }
    else if (parent_tag == "Modification")
    {
      actual_peptide_.mods.back().setMetaValue(name, data_value);
    }
    else if (parent_tag == "Compound")
    {
      actual_compound_.setMetaValue(name, data_value);
    }
    else if (parent_tag == "Protein")
    {
      actual_protein_.setMetaValue(name, data_value);
    }
    else if (parent_tag == "Configuration")
    {
      actual_configuration_.setMetaValue(name, data_value);
    }
    else if (parent_tag == "Prediction")
    {
      actual_prediction_.setMetaValue(name, data_value);
    }
    else if (parent_tag == "Interpretation")
    {
      actual_interpretation_.setMetaValue(name, data_value);
    }
    else if (parent_tag == "ValidationStatus")
    {
      actual_validation_.setMetaValue(name, data_value);
    }
    else if (parent_tag == "TargetList")
    {
      exp_->setTargetMetaValue(name, data_value);
    }
    else if (parent_tag == "Target")
    {
      actual_target_.setMetaValue(name, data_value);
    }
    else if (parent_tag == "Precursor")
    {
       if (parent_parent_tag == "Transition")
      {
        actual_transition_.setMetaValue(name, data_value);
      }
      if (parent_parent_tag == "Target")
      {
        actual_target_.setMetaValue(name, data_value);
      }
    }
    else if (parent_tag == "Product")
    {
      actual_transition_.setMetaValue(name, data_value);
    }
    else if (parent_tag == "SourceFile")
    {
      actual_sourcefile_.setMetaValue(name, data_value);
    } 
    else if (parent_tag == "Transition")
    {
      actual_transition_.setMetaValue(name, data_value);
    }
    else warning(LOAD, String("Unhandled userParam '") + name + "' in tag '" + parent_tag + "'.");
  }

  void TraMLHandler::writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent) const
  {
    std::vector<String> keys;
    meta.getKeys(keys);
    
    for (Size i = 0; i!=keys.size();++i)
    {
      os << String(2 * indent, ' ') << "<userParam name=\"" << keys[i] << "\" type=\"";
      
      DataValue d = meta.getMetaValue(keys[i]);
      //determine type
      if (d.valueType()==DataValue::INT_VALUE)
      {
        os << "xsd:integer";
      }
      else if (d.valueType()==DataValue::DOUBLE_VALUE)
      {
        os << "xsd:double";
      }
      else //string or lists are converted to string
      {
        os << "xsd:string";
      }
      os << "\" value=\"" << (String)(d) << "\"/>" << "\n";
    }
  }

  } //namespace Internal
} // namespace OpenMS


