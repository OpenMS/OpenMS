// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Andreas Bertsch, Hannes Roest $
// --------------------------------------------------------------------------

#include "OpenMS/CONCEPT/LogStream.h"
#include <OpenMS/FORMAT/HANDLERS/TraMLHandler.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <ostream>
#include <map>

namespace OpenMS::Internal
{


    TraMLHandler::TraMLHandler(const TargetedExperiment& exp, const String& filename, const String& version, const ProgressLogger& logger) :
      XMLHandler(filename, version),
      logger_(logger),
      exp_(nullptr),
      cexp_(&exp)
    {
      cv_.loadFromOBO("PI", File::find("/CV/psi-ms.obo"));
    }

    TraMLHandler::TraMLHandler(TargetedExperiment& exp, const String& filename, const String& version, const ProgressLogger& logger) :
      XMLHandler(filename, version),
      logger_(logger),
      exp_(&exp),
      cexp_(nullptr)
    {
      cv_.loadFromOBO("PI", File::find("/CV/psi-ms.obo"));
    }

    TraMLHandler::~TraMLHandler()
    = default;

    void TraMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      static const XMLCh* s_type = xercesc::XMLString::transcode("type");
      static const XMLCh* s_value = xercesc::XMLString::transcode("value");
      static const XMLCh* s_name = xercesc::XMLString::transcode("name");
      static const XMLCh* s_id = xercesc::XMLString::transcode("id");
      static const XMLCh* s_sequence = xercesc::XMLString::transcode("sequence");
      static const XMLCh* s_fullName = xercesc::XMLString::transcode("fullName");
      static const XMLCh* s_version = xercesc::XMLString::transcode("version");
      static const XMLCh* s_URI = xercesc::XMLString::transcode("URI");

      tag_ = sm_.convert(qname);
      open_tags_.push_back(tag_);

      static std::set<String> tags_to_ignore;
      if (tags_to_ignore.empty())
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
      if (open_tags_.size() > 1)
      {
        parent_tag = *(open_tags_.end() - 2);
      }
      String parent_parent_tag;
      if (open_tags_.size() > 2)
      {
        parent_parent_tag = *(open_tags_.end() - 3);
      }

      if (tag_ == "cvParam")
      {
        // These are here because of cppcheck
        static const XMLCh* s_accession = xercesc::XMLString::transcode("accession");
        static const XMLCh* s_unit_accession = xercesc::XMLString::transcode("unitAccession");
        static const XMLCh* s_unit_name = xercesc::XMLString::transcode("unitName");
        static const XMLCh* s_unit_cvref = xercesc::XMLString::transcode("unitCvRef");
        static const XMLCh* s_unit_ref = xercesc::XMLString::transcode("cvRef");

        String value, cv_ref, unit_accession, unit_name, unit_cv_ref;
        optionalAttributeAsString_(value, attributes, s_value);
        optionalAttributeAsString_(unit_accession, attributes, s_unit_accession);
        optionalAttributeAsString_(unit_name, attributes, s_unit_name);
        optionalAttributeAsString_(unit_cv_ref, attributes, s_unit_cvref);
        optionalAttributeAsString_(cv_ref, attributes, s_unit_ref);
        CVTerm::Unit unit(unit_accession, unit_name, unit_cv_ref);
        CVTerm cv_term(attributeAsString_(attributes, s_accession), 
                       attributeAsString_(attributes, s_name), cv_ref, value, unit);

        handleCVParam_(parent_parent_tag, parent_tag, cv_term);
        return;
      }
      else if (tag_ == "userParam")
      {
        String type = "";
        optionalAttributeAsString_(type, attributes, s_type);
        String value = "";
        optionalAttributeAsString_(value, attributes, s_value);
        handleUserParam_(parent_parent_tag, parent_tag, attributeAsString_(attributes, s_name), type, value);
      }
      else if (tag_ == "cv")
      {
        exp_->addCV(TargetedExperiment::CV(attributeAsString_(attributes, s_id), 
                                           attributeAsString_(attributes, s_fullName),
                                           attributeAsString_(attributes, s_version),
                                           attributeAsString_(attributes, s_URI)));
      }
      else if (tag_ == "Contact")
      {
        actual_contact_.id = attributeAsString_(attributes, s_id);
      }
      else if (tag_ == "Publication")
      {
        actual_publication_.id = attributeAsString_(attributes, s_id);
      }
      else if (tag_ == "Instrument")
      {
        actual_instrument_.id = attributeAsString_(attributes, s_id);
      }
      else if (tag_ == "Software")
      {
        actual_software_.setName(attributeAsString_(attributes, s_id));
        actual_software_.setVersion(attributeAsString_(attributes, s_version));
      }
      else if (tag_ == "Protein")
      {
        actual_protein_ = TargetedExperiment::Protein();
        actual_protein_.id = attributeAsString_(attributes, s_id);
      }
      else if (tag_ == "Peptide")
      {
        actual_peptide_ = TargetedExperiment::Peptide();
        actual_peptide_.id = attributeAsString_(attributes, s_id);
        actual_peptide_.sequence = attributeAsString_(attributes, s_sequence);
      }
      else if (tag_ == "Modification")
      {
        TargetedExperiment::Peptide::Modification mod;
        double avg_mass_delta(0), mono_mass_delta(0); // zero means no value
        optionalAttributeAsDouble_(avg_mass_delta, attributes, "averageMassDelta");
        optionalAttributeAsDouble_(mono_mass_delta, attributes, "monoisotopicMassDelta");
        mod.avg_mass_delta = avg_mass_delta;
        mod.mono_mass_delta = mono_mass_delta;

        mod.location = attributeAsInt_(attributes, "location") - 1; // TraML stores location starting with 1
        actual_peptide_.mods.push_back(mod);
      }
      else if (tag_ == "Compound")
      {
        actual_compound_ = TargetedExperiment::Compound();
        actual_compound_.id = attributeAsString_(attributes, s_id);
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
        if (optionalAttributeAsString_(id, attributes, s_id))
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
        actual_sourcefile_.setNativeIDType(attributeAsString_(attributes, s_id));
        actual_sourcefile_.setNameOfFile(attributeAsString_(attributes, s_name));
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
        if (optionalAttributeAsString_(id, attributes, s_id))
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
        actual_protein_.sequence = sm_.convert(chars);
        return;
      }
      return;
    }

    void TraMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      tag_ = sm_.convert(qname);

      //determine parent tag
      String parent_tag;
      if (open_tags_.size() > 1)
        parent_tag = *(open_tags_.end() - 2);
      String parent_parent_tag;
      if (open_tags_.size() > 2)
        parent_parent_tag = *(open_tags_.end() - 3);

      open_tags_.pop_back();

      static std::set<String> tags_to_ignore;
      if (tags_to_ignore.empty())
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
        actual_interpretation_ = TargetedExperiment::Interpretation();
      }
      else if (tag_ == "Prediction")
      {
        actual_transition_.setPrediction(actual_prediction_);
        actual_prediction_ = TargetedExperiment::Prediction();
      }
      else if (tag_ == "Configuration")
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
      logger_.startProgress(0, exp.getTransitions().size(), "storing TraML file");
      // int progress = 0;

      os << R"(<?xml version="1.0" encoding="UTF-8"?>)" << "\n";
      os << R"(<TraML version="1.0.0" xmlns="http://psi.hupo.org/ms/traml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/traml TraML1.0.0.xsd">)" << "\n";

      //--------------------------------------------------------------------------------------------
      // CV list
      //--------------------------------------------------------------------------------------------
      os << "  <cvList>" << "\n";

      if (exp.getCVs().empty())
      {
        os << R"(    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="unknown" URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo"/>)" << "\n"
           << R"(    <cv id="UO" fullName="Unit Ontology" version="unknown" URI="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo"/>)" << "\n";
      }
      else
      {
        for (std::vector<TargetedExperiment::CV>::const_iterator it = exp.getCVs().begin(); it != exp.getCVs().end(); ++it)
        {
          os << "    <cv id=\"" << writeXMLEscape(it->id) << "\" fullName=\"" << writeXMLEscape(it->fullname) << "\" version=\"" << writeXMLEscape(it->version) << "\" URI=\"" << writeXMLEscape(it->URI) << "\"/>" << "\n";
        }
      }
      os << "  </cvList>" << "\n";

      // source file list
      if (!exp.getSourceFiles().empty())
      {
        os << "  <SourceFileList>" << "\n";
        for (std::vector<SourceFile>::const_iterator it = exp.getSourceFiles().begin(); it != exp.getSourceFiles().end(); ++it)
        {
          os << "    <SourceFile id=\""
             << writeXMLEscape(it->getNativeIDType()) << "\" name=\""
             << writeXMLEscape(it->getNameOfFile()) << "\" location=\""
             << writeXMLEscape(it->getPathToFile()) << "\">"
             << "\n";
          writeCVParams_(os, *it, 3);
          writeUserParam_(os, (MetaInfoInterface) * it, 3);
          os << "    </SourceFile>" << "\n";
        }
        os << "  </SourceFileList>" << "\n";
      }

      // contact list
      if (!exp.getContacts().empty())
      {
        os << "  <ContactList>" << "\n";
        for (std::vector<TargetedExperiment::Contact>::const_iterator it = exp.getContacts().begin(); it != exp.getContacts().end(); ++it)
        {
          os << "    <Contact id=\"" << writeXMLEscape(it->id) << "\">" << "\n";
          writeCVParams_(os, *it, 3);
          writeUserParam_(os, (MetaInfoInterface) * it, 3);
          os << "    </Contact>" << "\n";
        }
        os << "  </ContactList>" << "\n";
      }

      // publication list
      if (!exp.getPublications().empty())
      {
        os << "  <PublicationList>"  << "\n";
        for (std::vector<TargetedExperiment::Publication>::const_iterator it = exp.getPublications().begin(); it != exp.getPublications().end(); ++it)
        {
          os << "    <Publication id=\"" << writeXMLEscape(it->id) << "\">" << "\n";
          writeCVParams_(os, *it, 3);
          writeUserParam_(os, (MetaInfoInterface) * it, 3);
          os << "    </Publication>" << "\n";
        }
        os << "  </PublicationList>" << "\n";
      }

      // instrument list
      if (!exp.getInstruments().empty())
      {
        os << "  <InstrumentList>" << "\n";
        for (std::vector<TargetedExperiment::Instrument>::const_iterator it = exp.getInstruments().begin(); it != exp.getInstruments().end(); ++it)
        {
          os << "    <Instrument id=\"" << writeXMLEscape(it->id) << "\">" << "\n";
          writeCVParams_(os, *it, 3);
          writeUserParam_(os, (MetaInfoInterface) * it, 3);
          os << "    </Instrument>" << "\n";
        }
        os << "  </InstrumentList>" << "\n";
      }

      // software list
      if (!exp.getSoftware().empty())
      {
        os << "  <SoftwareList>" << "\n";
        for (std::vector<Software>::const_iterator it = exp.getSoftware().begin(); it != exp.getSoftware().end(); ++it)
        {
          os << "    <Software id=\"" << writeXMLEscape(it->getName()) << "\" version=\"" << writeXMLEscape(it->getVersion()) << "\">" << "\n";
          writeCVParams_(os, *it, 3);
          writeUserParam_(os, (MetaInfoInterface) * it, 3);
          os << "    </Software>" << "\n";
        }
        os << "  </SoftwareList>" << "\n";
      }

      //--------------------------------------------------------------------------------------------
      // protein list
      //--------------------------------------------------------------------------------------------
      if (!exp.getProteins().empty())
      {
        os << "  <ProteinList>" << "\n";
        for (std::vector<TargetedExperiment::Protein>::const_iterator it = exp.getProteins().begin(); it != exp.getProteins().end(); ++it)
        {
          os << "    <Protein id=\"" << writeXMLEscape(it->id) << "\">" << "\n";
          writeCVParams_(os, *it, 3);
          writeUserParam_(os, (MetaInfoInterface) * it, 3);
          os << "      <Sequence>" << it->sequence << "</Sequence>" << "\n";
          os << "    </Protein>" << "\n";
        }
        os << "  </ProteinList>" << "\n";
      }

      //--------------------------------------------------------------------------------------------
      // compound list
      //--------------------------------------------------------------------------------------------
      ModificationsDB* mod_db = ModificationsDB::getInstance();
      if (exp.getCompounds().size()  + exp.getPeptides().size() > 0)
      {
        os << "  <CompoundList>" << "\n";
        std::vector<TargetedExperiment::Peptide> exp_peptides = exp.getPeptides();

        // 1. do peptides
        for (std::vector<TargetedExperiment::Peptide>::const_iterator it = exp_peptides.begin(); it != exp_peptides.end(); ++it)
        {
          os << "    <Peptide id=\"" << writeXMLEscape(it->id) << "\" sequence=\"" << it->sequence << "\">" << "\n";
          if (it->hasCharge())
          {
            os << R"(      <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value=")" <<  it->getChargeState() << "\"/>\n";
          }
          if (!it->getPeptideGroupLabel().empty())
          {
            os << R"(      <cvParam cvRef="MS" accession="MS:1000893" name="peptide group label" value=")" <<  it->getPeptideGroupLabel() << "\"/>\n";
          }
          if (it->getDriftTime() >= 0.0)
          {
            os << R"(      <cvParam cvRef="MS" accession="MS:1002476" name="ion mobility drift time" value=")" << it->getDriftTime() << "\" unitAccession=\"UO:0000028\" unitName=\"millisecond\" unitCvRef=\"UO\" />\n";
          }
          writeCVParams_(os,  *it, 3);
          writeUserParam_(os, (MetaInfoInterface) * it, 3);

          for (std::vector<String>::const_iterator rit = it->protein_refs.begin(); rit != it->protein_refs.end(); ++rit)
          {
            os << "      <ProteinRef ref=\"" << writeXMLEscape(*rit) << "\"/>" << "\n";
          }

          if (!it->mods.empty())
          {
            for (std::vector<TargetedExperiment::Peptide::Modification>::const_iterator 
                mit = it->mods.begin(); mit != it->mods.end(); ++mit)
            {
              os << "      <Modification";
              os << " location=\"" << mit->location + 1 << "\""; // TraML stores locations starting with 1
              if (mit->mono_mass_delta != 0)
              {
                os << " monoisotopicMassDelta=\"" << mit->mono_mass_delta << "\"";
              }
              if (mit->avg_mass_delta != 0)
              {
                os << " averageMassDelta=\"" << mit->avg_mass_delta << "\"";
              }
              os << ">\n";
              if (mit->unimod_id != -1)
              {
                // Get the name of the modifications from its unimod identifier (using getId)
                ResidueModification::TermSpecificity term_spec = ResidueModification::ANYWHERE;
                String residue = "";
                if (mit->location < 0)
                {
                  term_spec = ResidueModification::N_TERM;
                  if (!it->sequence.empty()) residue = it->sequence[0];
                }
                else if (Size(mit->location) >= it->sequence.size())
                {
                  term_spec = ResidueModification::C_TERM;
                  if (!it->sequence.empty()) residue = it->sequence[it->sequence.size() - 1];
                }
                else if (!it->sequence.empty())
                {
                  residue = it->sequence[mit->location];
                }
                const ResidueModification* rmod = mod_db->getModification("UniMod:" + String(mit->unimod_id), residue, term_spec);
                const String& modname = rmod->getId();
                os << R"(        <cvParam cvRef="UNIMOD" accession="UNIMOD:)" << mit->unimod_id
                  << "\" name=\"" << modname << "\"/>\n";
              }

              writeCVParams_(os, *mit, 4);
              writeUserParam_(os, (MetaInfoInterface) * mit, 4);
              os << "      </Modification>\n";
            }
          }

          if (!it->rts.empty())
          {
            os << "      <RetentionTimeList>\n";
            for (std::vector<TargetedExperiment::RetentionTime>::const_iterator rit = it->rts.begin(); rit != it->rts.end(); ++rit)
            {
              writeRetentionTime_(os, *rit);
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

        // 2. do compounds
        for (std::vector<TargetedExperiment::Compound>::const_iterator it = exp.getCompounds().begin(); it != exp.getCompounds().end(); ++it)
        {
          os << "    <Compound id=\"" << writeXMLEscape(it->id) << "\">" << "\n";

          if (it->hasCharge())
          {
            os << R"(      <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value=")" <<  it->getChargeState() << "\"/>\n";
          }
          if (it->theoretical_mass > 0.0)
          {
            os << R"(      <cvParam cvRef="MS" accession="MS:1001117" name="theoretical mass" value=")" << 
            it->theoretical_mass << "\" unitCvRef=\"UO\" unitAccession=\"UO:0000221\" unitName=\"dalton\"/>\n";
          }
          if (!it->molecular_formula.empty())
          {
            os << R"(      <cvParam cvRef="MS" accession="MS:1000866" name="molecular formula" value=")" << 
            it->molecular_formula << "\"/>\n";
          }
          if (!it->smiles_string.empty())
          {
            os << R"(      <cvParam cvRef="MS" accession="MS:1000868" name="SMILES string" value=")" << 
            it->smiles_string << "\"/>\n";
          }

          writeCVParams_(os, *it, 3);
          writeUserParam_(os, (MetaInfoInterface) * it, 3);

          if (!it->rts.empty())
          {
            os << "      <RetentionTimeList>\n";
            for (std::vector<TargetedExperiment::RetentionTime>::const_iterator rit = it->rts.begin(); rit != it->rts.end(); ++rit)
            {
              writeRetentionTime_(os, *rit);
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
      if (!exp.getTransitions().empty())
      {
        int progress = 0;

        os << "  <TransitionList>" << "\n";
        for (std::vector<ReactionMonitoringTransition>::const_iterator it = exp.getTransitions().begin(); it != exp.getTransitions().end(); ++it)
        {
          logger_.setProgress(progress++);
          os << "    <Transition";
          os << " id=\"" << writeXMLEscape(it->getName()) << "\"";

          if (!it->getPeptideRef().empty())
          {
            os << " peptideRef=\"" << writeXMLEscape(it->getPeptideRef()) << "\"";
          }

          if (!it->getCompoundRef().empty())
          {
            os << " compoundRef=\"" << writeXMLEscape(it->getCompoundRef()) << "\"";
          }
          os << ">" << "\n";

          // Precursor occurs exactly once (is required according to schema).
          // CV term MS:1000827 MUST be supplied for the TransitionList path
          os << "      <Precursor>" << "\n";
          os << R"(        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value=")" <<
            precisionWrapper(it->getPrecursorMZ()) << "\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n";
          if (it->hasPrecursorCVTerms())
          {
            writeCVParams_(os, it->getPrecursorCVTermList(), 4);
            writeUserParam_(os, (MetaInfoInterface)it->getPrecursorCVTermList(), 4);
          }
          os << "      </Precursor>" << "\n";

          for (ProductListType::const_iterator prod_it = it->getIntermediateProducts().begin();
               prod_it != it->getIntermediateProducts().end(); ++prod_it)
          {
            os << "      <IntermediateProduct>" << "\n";
            writeProduct_(os, prod_it);
            os << "      </IntermediateProduct>" << "\n";
          }

          // Product is required
          os << "      <Product>" << "\n";
          ProductListType dummy_vect;
          dummy_vect.push_back(it->getProduct());
          writeProduct_(os, dummy_vect.begin());
          os << "      </Product>" << "\n";

          const TargetedExperimentHelper::RetentionTime rit = it->getRetentionTime();
          if (!rit.getCVTerms().empty())
          {
            writeRetentionTime_(os, rit);
          }

          if (it->hasPrediction())
          {
            os << "      <Prediction softwareRef=\"" << writeXMLEscape(it->getPrediction().software_ref) << "\"";
            if (!it->getPrediction().contact_ref.empty())
            {
              os << " contactRef=\"" << writeXMLEscape(it->getPrediction().contact_ref) << "\"";
            }
            os << ">" << "\n";
            writeCVParams_(os, it->getPrediction(), 4);
            writeUserParam_(os, (MetaInfoInterface)it->getPrediction(), 4);
            os << "      </Prediction>" << "\n";
          }

          writeCVParams_(os, *it, 3);
          // Special CV Params
          if (it->getLibraryIntensity() > -100)
          {
            os << R"(      <cvParam cvRef="MS" accession="MS:1001226" name="product ion intensity" value=")" <<  it->getLibraryIntensity() << "\"/>\n";
          }
          if (it->getDecoyTransitionType() != ReactionMonitoringTransition::UNKNOWN)
          {
            if (it->getDecoyTransitionType() == ReactionMonitoringTransition::TARGET)
            {
              os << "      <cvParam cvRef=\"MS\" accession=\"MS:1002007\" name=\"target SRM transition\"/>\n";
            }
            else if (it->getDecoyTransitionType() == ReactionMonitoringTransition::DECOY)
            {
              os << "      <cvParam cvRef=\"MS\" accession=\"MS:1002008\" name=\"decoy SRM transition\"/>\n";
            }
          }

          // Output transition type (only write if non-default, otherwise assume default)
          // Default is: true, false, true
          // NOTE: do not change that, the same default is implicitly assumed in ReactionMonitoringTransition
          if (!it->isDetectingTransition())
          {
            os << "      <userParam name=\"detecting_transition\" type=\"xsd:boolean\" value=\"false\"/>\n";
          }
          if (it->isIdentifyingTransition())
          {
            os << "      <userParam name=\"identifying_transition\" type=\"xsd:boolean\" value=\"true\"/>\n";
          }
          if (!it->isQuantifyingTransition())
          {
            os << "      <userParam name=\"quantifying_transition\" type=\"xsd:boolean\" value=\"false\"/>\n";
          }

          writeUserParam_(os, (MetaInfoInterface) * it, 3);

          os << "    </Transition>" << "\n";
        }
        os << "  </TransitionList>" << "\n";

      }

      if (!exp.getIncludeTargets().empty() || !exp.getExcludeTargets().empty())
      {
        os << "  <TargetList>" << "\n";
        writeCVParams_(os, exp.getTargetCVTerms(), 2);
        writeUserParam_(os, (MetaInfoInterface)exp.getTargetCVTerms(), 2);

        if (!exp.getIncludeTargets().empty())
        {
          os << "    <TargetIncludeList>" << "\n";
          for (std::vector<IncludeExcludeTarget>::const_iterator it = exp.getIncludeTargets().begin(); it != exp.getIncludeTargets().end(); ++it)
          {
            writeTarget_(os, it);
          }
          os << "    </TargetIncludeList>" << "\n";
        }

        if (!exp.getExcludeTargets().empty())
        {
          os << "    <TargetExcludeList>" << "\n";
          for (std::vector<IncludeExcludeTarget>::const_iterator it = exp.getExcludeTargets().begin(); it != exp.getExcludeTargets().end(); ++it)
          {
            writeTarget_(os, it);
          }
          os << "    </TargetExcludeList>" << "\n";
        }

        os << "  </TargetList>" << "\n";
      }

      os << "</TraML>" << "\n";
      logger_.endProgress();
      return;
    }

    void TraMLHandler::writeRetentionTime_(std::ostream& os, const TargetedExperimentHelper::RetentionTime& rt) const
    {
      const TargetedExperimentHelper::RetentionTime* rit = &rt;
      os << "        <RetentionTime";
      if (!rit->software_ref.empty())
      {
        os << " softwareRef=\"" << writeXMLEscape(rit->software_ref) << "\"";
      }
      os << ">" << "\n";

      if (rit->isRTset())
      {
        if (rit->retention_time_type == TargetedExperimentHelper::RetentionTime::RTType::LOCAL)
        {
          os << R"(          <cvParam cvRef="MS" accession="MS:1000895" name="local retention time" value=")" << rit->getRT() << "\"";
        }
        else if (rit->retention_time_type == TargetedExperimentHelper::RetentionTime::RTType::NORMALIZED)
        {
          os << R"(          <cvParam cvRef="MS" accession="MS:1000896" name="normalized retention time" value=")" << rit->getRT() << "\"";
        }
        else if (rit->retention_time_type == TargetedExperimentHelper::RetentionTime::RTType::PREDICTED)
        {
          os << R"(          <cvParam cvRef="MS" accession="MS:1000897" name="predicted retention time" value=")" << rit->getRT() << "\"";
        }
        else if (rit->retention_time_type == TargetedExperimentHelper::RetentionTime::RTType::HPINS)
        {
          os << R"(          <cvParam cvRef="MS" accession="MS:1000902" name="H-PINS retention time normalization standard" value=")" << rit->getRT() << "\"";
        }
        else if (rit->retention_time_type == TargetedExperimentHelper::RetentionTime::RTType::IRT)
        {
          os << R"(          <cvParam cvRef="MS" accession="MS:1002005" name="iRT retention time normalization standard" value=")" << rit->getRT() << "\"";
        }
        else
        {
          os << R"(          <cvParam cvRef="MS" accession="MS:1000895" name="local retention time" value=")" << rit->getRT() << "\"";
        }
      }

      // write units (minute, second or none)
      if ( rit->retention_time_unit == TargetedExperimentHelper::RetentionTime::RTUnit::SECOND) //seconds
      {
        os << " unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"/>\n";
      }
      else if ( rit->retention_time_unit == TargetedExperimentHelper::RetentionTime::RTUnit::MINUTE) //minutes
      {
        os << " unitCvRef=\"UO\" unitAccession=\"UO:0000031\" unitName=\"minute\"/>\n";
      }
      else
      {
        os << "/>\n";
      }

      writeCVParams_(os, *rit, 5);
      writeUserParam_(os, (MetaInfoInterface) * rit, 5);
      os << "        </RetentionTime>" << "\n";
    }

    void TraMLHandler::writeTarget_(std::ostream& os, const std::vector<IncludeExcludeTarget>::const_iterator& it) const
    {
      os << "      <Target id=\"" << writeXMLEscape(it->getName()) << "\"";
      if (!it->getPeptideRef().empty())
      {
        os << " peptideRef=\"" << writeXMLEscape(it->getPeptideRef()) << "\"";
      }
      if (!it->getCompoundRef().empty())
      {
        os << " compoundRef=\"" << writeXMLEscape(it->getCompoundRef()) << "\"";
      }
      os << ">\n";
      os << "        <Precursor>\n";
      writeCVParams_(os, it->getPrecursorCVTermList(), 5);
      writeUserParam_(os, (MetaInfoInterface)it->getPrecursorCVTermList(), 5);
      os << "        </Precursor>\n";

      const IncludeExcludeTarget::RetentionTime* rit = &it->getRetentionTime();
      if (!rit->getCVTerms().empty())
      {
        writeRetentionTime_(os, *rit);
      }

      if (!it->getConfigurations().empty())
      {
        os << "        <ConfigurationList>\n";
        for (auto config_it = it->getConfigurations().begin(); config_it != it->getConfigurations().end(); ++config_it)
        {
          writeConfiguration_(os, config_it);
        }
        os << "        </ConfigurationList>\n";
      }

      // TODO  : add cv/userparams for Target
      os << "      </Target>" << "\n";

    }

    void TraMLHandler::writeProduct_(std::ostream& os, const std::vector<ReactionMonitoringTransition::Product>::const_iterator& prod_it) const
    {
      if (prod_it->hasCharge())
      {
        os << R"(        <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value=")" <<  prod_it->getChargeState() << "\"/>\n";
      }
      if (prod_it->getMZ() > 0)
      {
        os << R"(        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value=")" <<  
          prod_it->getMZ() << "\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n";
      }
      writeCVParams_(os, *prod_it, 4);
      writeUserParam_(os, (MetaInfoInterface) * prod_it, 4);

      if (!prod_it->getInterpretationList().empty())
      {
        os << "        <InterpretationList>" << "\n";
        for (std::vector<TargetedExperiment::Interpretation>::const_iterator inter_it = prod_it->getInterpretationList().begin(); 
            inter_it != prod_it->getInterpretationList().end(); ++inter_it)
        {
          os << "          <Interpretation>" << "\n";
          if (inter_it->ordinal > 0)
          {
            os << R"(            <cvParam cvRef="MS" accession="MS:1000903" name="product ion series ordinal" value=")" << 
              (int)inter_it->ordinal << "\"/>\n";
          }
          if (inter_it->rank > 0)
          {
            os << R"(            <cvParam cvRef="MS" accession="MS:1000926" name="product interpretation rank" value=")" << 
              (int)inter_it->rank << "\"/>\n";
          }

          // Ion Type
          switch (inter_it->iontype)
          {
            case Residue::AIon:
              os << "            <cvParam cvRef=\"MS\" accession=\"MS:1001229\" name=\"frag: a ion\"/>\n";
              break;
            case Residue::BIon:
              os << "            <cvParam cvRef=\"MS\" accession=\"MS:1001224\" name=\"frag: b ion\"/>\n";
              break;
            case Residue::CIon:
              os << "            <cvParam cvRef=\"MS\" accession=\"MS:1001231\" name=\"frag: c ion\"/>\n";
              break;
            case Residue::XIon:
              os << "            <cvParam cvRef=\"MS\" accession=\"MS:1001228\" name=\"frag: x ion\"/>\n";
              break;
            case Residue::YIon:
              os << "            <cvParam cvRef=\"MS\" accession=\"MS:1001220\" name=\"frag: y ion\"/>\n";
              break;
            case Residue::ZIon:
              os << "            <cvParam cvRef=\"MS\" accession=\"MS:1001230\" name=\"frag: z ion\"/>\n";
              break;
            case Residue::Precursor:
              os << "            <cvParam cvRef=\"MS\" accession=\"MS:1001523\" name=\"frag: precursor ion\"/>\n";
              break;
            case Residue::BIonMinusH20:
              os << "            <cvParam cvRef=\"MS\" accession=\"MS:1001222\" name=\"frag: b ion - H2O\"/>\n";
              break;
            case Residue::YIonMinusH20:
              os << "            <cvParam cvRef=\"MS\" accession=\"MS:1001223\" name=\"frag: y ion - H2O\"/>\n";
              break;
            case Residue::BIonMinusNH3:
              os << "            <cvParam cvRef=\"MS\" accession=\"MS:1001232\" name=\"frag: b ion - NH3\"/>\n";
              break;
            case Residue::YIonMinusNH3:
              os << "            <cvParam cvRef=\"MS\" accession=\"MS:1001233\" name=\"frag: y ion - NH3\"/>\n";
              break;
            case Residue::NonIdentified:
              os << "            <cvParam cvRef=\"MS\" accession=\"MS:1001240\" name=\"non-identified ion\"/>\n";
              break;
            case Residue::Unannotated:
              // means no annotation and no input cvParam - to write out a cvParam, use Residue::NonIdentified
              break;
            // invalid values
            case Residue::Zp1Ion: OPENMS_LOG_ERROR << "Zp1 ions not supported. Ignoring." << std::endl; break;
            case Residue::Zp2Ion: OPENMS_LOG_ERROR << "Zp2 ions not supported. Ignoring." << std::endl; break;
            case Residue::Full: break;
            case Residue::Internal: break;
            case Residue::NTerminal: break;
            case Residue::CTerminal: break;
            case Residue::SizeOfResidueType:
              break;
          }
          writeCVParams_(os, *inter_it, 6);
          writeUserParam_(os, (MetaInfoInterface) * inter_it, 6);
          os << "          </Interpretation>" << "\n";
        }
        os << "        </InterpretationList>" << "\n";
      }
      if (!prod_it->getConfigurationList().empty())
      {
        os << "        <ConfigurationList>" << "\n";
        for (ConfigurationListType::const_iterator config_it = prod_it->getConfigurationList().begin(); config_it != prod_it->getConfigurationList().end(); ++config_it)
        {
          writeConfiguration_(os, config_it);
        }
        os << "        </ConfigurationList>" << "\n";
      }
    }

    void TraMLHandler::writeConfiguration_(std::ostream& os, const std::vector<ReactionMonitoringTransition::Configuration>::const_iterator& cit) const
    {
      os << "          <Configuration instrumentRef=\"" << writeXMLEscape(cit->instrument_ref) << "\"";
      if (!cit->contact_ref.empty())
      {
        os << " contactRef=\"" << writeXMLEscape(cit->contact_ref) << "\"";
      }
      os << ">" << "\n";

      writeCVParams_(os, *cit, 6);
      writeUserParam_(os, (MetaInfoInterface) * cit, 6);
      if (!cit->validations.empty())
      {
        for (std::vector<CVTermList>::const_iterator iit = cit->validations.begin(); iit != cit->validations.end(); ++iit)
        {
          if (!iit->empty())
          {
            os << "            <ValidationStatus>" << "\n";
            writeCVParams_(os, *iit, 7);
            writeUserParam_(os, (MetaInfoInterface) * iit, 7);
            os << "            </ValidationStatus>" << "\n";
          }
        }
      }

      os << "          </Configuration>" << "\n";
    }

    void TraMLHandler::handleCVParam_(const String& parent_parent_tag, const String& parent_tag, const CVTerm& cv_term)
    {
      //Error checks of CV values
      const String& accession = cv_term.getAccession();
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
        if (parsed_name != correct_name)
        {
          warning(LOAD, String("Name of CV term not correct: '") + term.id + " - " + parsed_name + "' should be '" + correct_name + "'");
        }
        if (term.obsolete)
        {
          warning(LOAD, String("Obsolete CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "'.");
          //values used in wrong places and wrong value types
          String value = cv_term.getValue().toString();
          if (!value.empty())
          {
            if (term.xref_type == ControlledVocabulary::CVTerm::NONE)
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
                catch (Exception::ConversionError&)
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
                catch (Exception::ConversionError&)
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
                catch (Exception::ParseError&)
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
          else if (term.xref_type != ControlledVocabulary::CVTerm::NONE && term.xref_type != ControlledVocabulary::CVTerm::XSD_STRING)
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
        // Note: we have to be prepared to have multiple CV terms for the same
        // RT, some indicating the unit, some indicating the type of RT

        // MAY supply a *child* term of MS:1000915 (retention time window attribute) one or more times
        //   e.g.: MS:1000916 (retention time window lower offset)
        //   e.g.: MS:1000917 (retention time window upper offset)
        //   e.g.: MS:1001907 (retention time window width)
        // MAY supply a *child* term of MS:1000901 (retention time normalization standard) only once
        //   e.g.: MS:1000902 (H-PINS retention time normalization standard)
        //   e.g.: MS:1002005 (iRT retention time normalization standard)
        // MAY supply a *child* term of MS:1000894 (retention time) one or more times
        //   e.g.: MS:1000895 (local retention time)
        //   e.g.: MS:1000896 (normalized retention time)
        //   e.g.: MS:1000897 (predicted retention time)

        if ( cv_term.getUnit().accession == "UO:0000010") //seconds
        {
          actual_rt_.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
        }
        else if ( cv_term.getUnit().accession == "UO:0000031") //minutes
        {
          actual_rt_.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::MINUTE;
        }
        else if (actual_rt_.retention_time_unit == TargetedExperimentHelper::RetentionTime::RTUnit::SIZE_OF_RTUNIT) // do not overwrite previous data
        {
          actual_rt_.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::UNKNOWN;
        }

        if (cv_term.getAccession() == "MS:1000895") // local RT
        {
          actual_rt_.setRT(cv_term.getValue().toString().toDouble());
          actual_rt_.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::LOCAL;
        }
        else if (cv_term.getAccession() == "MS:1000896") // normalized RT
        {
          actual_rt_.setRT(cv_term.getValue().toString().toDouble());
          actual_rt_.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::NORMALIZED;
        }
        else if (cv_term.getAccession() == "MS:1000897") // predicted RT
        {
          actual_rt_.setRT(cv_term.getValue().toString().toDouble());
          actual_rt_.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::PREDICTED;
        }
        else if (cv_term.getAccession() == "MS:1000902") // H-PINS
        {
          if (!cv_term.getValue().toString().empty()) actual_rt_.setRT(cv_term.getValue().toString().toDouble());
          actual_rt_.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::HPINS;
        }
        else if (cv_term.getAccession() == "MS:1002005") // iRT
        {
          if (!cv_term.getValue().toString().empty()) actual_rt_.setRT(cv_term.getValue().toString().toDouble());
          actual_rt_.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::IRT;
        }
        // else if (cv_term.getAccession() == "MS:1000916") // RT lower offset
        // {
        //   actual_rt_.retention_time_lower = cv_term.getValue().toString().toDouble();
        // }
        // else if (cv_term.getAccession() == "MS:1000917") // RT upper offset
        // {
        //   actual_rt_.retention_time_upper = cv_term.getValue().toString().toDouble();
        // }
        // else if (cv_term.getAccession() == "MS:1001907") // RT window width
        // {
        //   actual_rt_.retention_time_width = cv_term.getValue().toString().toDouble();
        // }
        else
        {
          warning(LOAD, String("The CV term '" + cv_term.getAccession() + "' - '" +
                cv_term.getName() + "' used in tag '" + parent_tag + "' is currently not supported!"));
          actual_rt_.addCVTerm(cv_term);
        }
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
        else if (cv_term.getAccession() == "MS:1002476")
        {
          actual_peptide_.setDriftTime(cv_term.getValue().toString().toDouble());
        }
        else
        {
          actual_peptide_.addCVTerm(cv_term);
        }
      }
      else if (parent_tag == "Modification")
      {
        // if we find a CV term that starts with UniMod, chances are we can use
        // the UniMod accession number to identify the modification
        if (cv_term.getAccession().size() > 7 && cv_term.getAccession().prefix(7).toLower() == String("unimod:"))
        {
          // check for Exception::ConversionError ?
          actual_peptide_.mods.back().unimod_id = cv_term.getAccession().substr(7).toInt();
        }
        else
        {
          actual_peptide_.mods.back().addCVTerm(cv_term);
        }

      }
      else if (parent_tag == "Compound")
      {
        if (cv_term.getAccession() == "MS:1001117")
        {
          actual_compound_.theoretical_mass = cv_term.getValue().toString().toDouble();
        }
        else if (cv_term.getAccession() == "MS:1000866")
        {
          actual_compound_.molecular_formula = cv_term.getValue().toString();
        }
        else if (cv_term.getAccession() == "MS:1000868")
        {
          actual_compound_.smiles_string = cv_term.getValue().toString();
        }
        else if (cv_term.getAccession() == "MS:1000041")
        {
          actual_compound_.setChargeState(cv_term.getValue().toString().toInt());
        }
        else if (cv_term.getAccession() == "MS:1002476")
        {
          actual_peptide_.setDriftTime(cv_term.getValue().toString().toDouble());
        }
        else
        {
          actual_compound_.addCVTerm(cv_term);
        }
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

        ////
        ////    enum ResidueType
        ////    {
        ////      Full = 0,       // with N-terminus and C-terminus
        ////      Internal,       // internal, without any termini
        ////      NTerminal,      // only N-terminus
        ////      CTerminal,      // only C-terminus
        ////      AIon,           // MS:1001229 N-terminus up to the C-alpha/carbonyl carbon bond
        ////      BIon,           // MS:1001224 N-terminus up to the peptide bond
        ////      CIon,           // MS:1001231 N-terminus up to the amide/C-alpha bond
        ////      XIon,           // MS:1001228 amide/C-alpha bond up to the C-terminus
        ////      YIon,           // MS:1001220 peptide bond up to the C-terminus
        ////      ZIon,           // MS:1001230 C-alpha/carbonyl carbon bond
        ////      Precursor,      // MS:1001523 Precursor ion
        ////      BIonMinusH20,   // MS:1001222 b ion without water
        ////      YIonMinusH20,   // MS:1001223 y ion without water
        ////      BIonMinusNH3,   // MS:1001232 b ion without ammonia
        ////      YIonMinusNH3,   // MS:1001233 y ion without ammonia
        ////      Unannotated,    // unknown annotation
        ////      SizeOfResidueType
        ////    };

        if (cv_term.getAccession() == "MS:1000903")
        {
          // name: product ion series ordinal
          // def: "The ordinal of the fragment within a specified ion series. (e.g. 8 for a y8 ion)." [PSI:PI]
          actual_interpretation_.ordinal = cv_term.getValue().toString().toInt();
        }
        else if (cv_term.getAccession() == "MS:1000926")
        {
          // name: product interpretation rank
          // def: "The integer rank given an interpretation of an observed product ion. For example, if y8 is selected as the most likely interpretation of a peak, then it is assigned a rank of 1." [PSI:MS]
          actual_interpretation_.rank = cv_term.getValue().toString().toInt();
        }
        else if (cv_term.getAccession() == "MS:1001229")
        {
          actual_interpretation_.iontype = TargetedExperiment::IonType::AIon;
        }
        else if (cv_term.getAccession() == "MS:1001224")
        {
          actual_interpretation_.iontype = TargetedExperiment::IonType::BIon;
        }
        else if (cv_term.getAccession() == "MS:1001231")
        {
          actual_interpretation_.iontype = TargetedExperiment::IonType::CIon;
        }
        else if (cv_term.getAccession() == "MS:1001228")
        {
          actual_interpretation_.iontype = TargetedExperiment::IonType::XIon;
        }
        else if (cv_term.getAccession() == "MS:1001220")
        {
          actual_interpretation_.iontype = TargetedExperiment::IonType::YIon;
        }
        else if (cv_term.getAccession() == "MS:1001230")
        {
          actual_interpretation_.iontype = TargetedExperiment::IonType::ZIon;
        }
        else if (cv_term.getAccession() == "MS:1001523")
        {
          actual_interpretation_.iontype = TargetedExperiment::IonType::Precursor;
        }
        else if (cv_term.getAccession() == "MS:1001222")
        {
          actual_interpretation_.iontype = TargetedExperiment::IonType::BIonMinusH20;
        }
        else if (cv_term.getAccession() == "MS:1001223")
        {
          actual_interpretation_.iontype = TargetedExperiment::IonType::YIonMinusH20;
        }
        else if (cv_term.getAccession() == "MS:1001232")
        {
          actual_interpretation_.iontype = TargetedExperiment::IonType::BIonMinusNH3;
        }
        else if (cv_term.getAccession() == "MS:1001233")
        {
          actual_interpretation_.iontype = TargetedExperiment::IonType::YIonMinusNH3;
        }
        else if (cv_term.getAccession() == "MS:1001240")
        {
          actual_interpretation_.iontype = TargetedExperiment::IonType::NonIdentified;
        }
        else
        {
          actual_interpretation_.addCVTerm(cv_term);
        }
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
        if (cv_term.getAccession() == "MS:1000041")
        {
          actual_product_.setChargeState(cv_term.getValue().toString().toDouble());
        }
        else if (cv_term.getAccession() == "MS:1000827")
        {
          actual_product_.setMZ(cv_term.getValue().toString().toDouble());
        }
        else
        {
          actual_product_.addCVTerm(cv_term);
        }
      }
      else if (parent_tag == "Product")
      {
        if (cv_term.getAccession() == "MS:1000041")
        {
          actual_product_.setChargeState(cv_term.getValue().toString().toDouble());
        }
        else if (cv_term.getAccession() == "MS:1000827")
        {
          actual_product_.setMZ(cv_term.getValue().toString().toDouble());
        }
        else
        {
          actual_product_.addCVTerm(cv_term);
        }
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
        warning(LOAD, String("The CV term '" + cv_term.getAccession() + "' - '" +
              cv_term.getName() + "' used in tag '" + parent_tag + "' could not be handled, ignoring it!"));
      }
      return;
    }

    void TraMLHandler::handleUserParam_(const String& parent_parent_tag, const String& parent_tag, const String& name, const String& type, const String& value)
    {
      // create a DataValue that contains the data in the right type
      DataValue data_value = fromXSDString(type, value);

      // find the right MetaInfoInterface
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
        // see xsd:boolean reference (http://books.xmlschemata.org/relaxng/ch19-77025.html)
        // The value space of xsd:boolean is true and false. Its lexical space
        // accepts true, false, and also 1 (for true) and 0 (for false).
        if (name == "detecting_transition")
        {
          actual_transition_.setDetectingTransition((value == "true" || value == "1"));
        }
        else if (name == "identifying_transition")
        {
          actual_transition_.setIdentifyingTransition((value == "true" || value == "1"));
        }
        else if (name == "quantifying_transition")
        {
          actual_transition_.setQuantifyingTransition((value == "true" || value == "1"));
        }
        else
        {
          actual_transition_.setMetaValue(name, data_value);
        }
      }
      else
      {
        warning(LOAD, String("Unhandled userParam '") + name + "' in tag '" + parent_tag + "'.");
      }
    }

    void TraMLHandler::writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent) const
    {
      std::vector<String> keys;
      meta.getKeys(keys);

      for (Size i = 0; i != keys.size(); ++i)
      {
        os << String(2 * indent, ' ') << "<userParam name=\"" << writeXMLEscape(keys[i]) << "\" type=\"";

        const DataValue& d = meta.getMetaValue(keys[i]);
        //determine type
        if (d.valueType() == DataValue::INT_VALUE)
        {
          os << "xsd:integer";
        }
        else if (d.valueType() == DataValue::DOUBLE_VALUE)
        {
          os << "xsd:double";
        }
        else //string or lists are converted to string
        {
          os << "xsd:string";
        }
        os << "\" value=\"" << writeXMLEscape((String)(d)) << "\"/>" << "\n";
      }
    }

    void TraMLHandler::writeCVParams_(std::ostream & os, const CVTermList & cv_terms, UInt indent) const
    {
      writeCVList_(os, cv_terms.getCVTerms(), indent);
    }

    void TraMLHandler::writeCVParams_(std::ostream & os, const CVTermListInterface & cv_terms, UInt indent) const
    {
      writeCVList_(os, cv_terms.getCVTerms(), indent);
    }

    void TraMLHandler::writeCVList_(std::ostream & os, const std::map<String, std::vector<CVTerm>> & cv_terms, UInt indent) const
    {
      for (std::map<String, std::vector<CVTerm> >::const_iterator it = cv_terms.begin();
           it != cv_terms.end(); ++it)
      {
        for (const CVTerm& cit : it->second)
        {
          os << String(2 * indent, ' ') << "<cvParam cvRef=\"" << cit.getCVIdentifierRef() << "\" accession=\"" << cit.getAccession() << "\" name=\"" << cit.getName() << "\"";
          if (cit.hasValue() && !cit.getValue().isEmpty() && !cit.getValue().toString().empty())
          {
            os << " value=\"" << cit.getValue().toString() << "\"";
          }

          if (cit.hasUnit())
          {
            os << " unitCvRef=\"" << cit.getUnit().cv_ref << "\" unitAccession=\"" << cit.getUnit().accession << "\" unitName=\"" << cit.getUnit().name << "\"";
          }
          os << "/>" << "\n";
        }
      }
    }
} // namespace OpenMS //namespace Internal

