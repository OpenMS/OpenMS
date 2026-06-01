// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/ProtXMLFile.h>

using namespace std;

namespace OpenMS
{

  ProtXMLFile::ProtXMLFile() :
    XMLHandler("", "1.2"),
    XMLFile("/SCHEMAS/protXML_v6.xsd", "6.0")
  {
  }

  void ProtXMLFile::load(const String& filename, ProteinIdentification& protein_ids, PeptideIdentification& peptide_ids)
  {
    //Filename for error messages in XMLHandler
    file_ = filename;

    resetMembers_();

    // reset incoming data
    protein_ids = ProteinIdentification();
    peptide_ids = PeptideIdentification();

    // remember data link while parsing
    prot_id_ = &protein_ids;
    pep_id_ = &peptide_ids;

    parse_(filename, this);
  }

  void ProtXMLFile::store(const String& /*filename*/, const ProteinIdentification& /*protein_ids*/, const PeptideIdentification& /*peptide_ids*/, const String& /*document_id*/)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    // resetMembers_();
  }

  /// reset members
  void ProtXMLFile::resetMembers_()
  {
    prot_id_ = nullptr;
    pep_id_ = nullptr;
    pep_hit_ = nullptr;
    protein_group_ = ProteinGroup();
  }

  void ProtXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
  {
    String tag = sm_.convert(qname);

    if (tag == "protein_summary_header")
    {
      String db = attributeAsString_(attributes, "reference_database");
      String enzyme = attributeAsString_(attributes, "sample_enzyme");
      ProteinIdentification::SearchParameters sp = prot_id_->getSearchParameters();
      sp.db = db;
      // find a matching enzyme name
      sp.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme));
      prot_id_->setSearchParameters(sp);
      prot_id_->setScoreType("ProteinProphet probability");
      prot_id_->setHigherScoreBetter(true);
      pep_id_->setScoreType("ProteinProphet probability");
      pep_id_->setHigherScoreBetter(true);
    }
    // identifier for Protein & PeptideIdentification
    // <program_details analysis="proteinprophet" time="2009-11-29T18:30:03" ...
    if (tag == "program_details")
    {
      String analysis = attributeAsString_(attributes, "analysis");
      String time = attributeAsString_(attributes, "time");
      String version = attributeAsString_(attributes, "version");

      DateTime date;
      date.set(time);

      if (!date.isValid())
      {
        OPENMS_LOG_WARN << "Warning: Cannot parse 'time'='" << time << "'.\n";
      }
      prot_id_->setDateTime(date);
      prot_id_->setSearchEngine(analysis);
      prot_id_->setSearchEngineVersion(version);
      String id = String(UniqueIdGenerator::getUniqueId()); // was: analysis + "_" + time;
      prot_id_->setIdentifier(id);
      pep_id_->setIdentifier(id);
    }

    if (tag == "protein_group")
    {
      // we group all <protein>'s and <indistinguishable_protein>'s in our
      // internal group structure
      protein_group_ = ProteinGroup();
      protein_group_.probability = attributeAsDouble_(attributes, "probability");
    }
    else if (tag == "protein")
    {
      // usually there will be just one <protein> per <protein_group>, but more
      // are possible; each <protein> is distinguishable from the other, we
      // nevertheless group them

      String protein_name = attributeAsString_(attributes, "protein_name");
      // open new "indistinguishable" group:
      prot_id_->insertIndistinguishableProteins(ProteinGroup());
      registerProtein_(protein_name); // create new protein

      // fill protein with life
      double pc_coverage;
      if (optionalAttributeAsDouble_(pc_coverage, attributes, "percent_coverage"))
      {
        prot_id_->getHits().back().setCoverage(pc_coverage);
      }
      else
      {
        OPENMS_LOG_WARN << "Required attribute 'percent_coverage' missing\n";
      }
      prot_id_->getHits().back().setScore(attributeAsDouble_(attributes, "probability"));

    }
    else if (tag == "indistinguishable_protein")
    {
      String protein_name = attributeAsString_(attributes, "protein_name");
      // current last protein is from the same "indistinguishable" group:
      double score = prot_id_->getHits().back().getScore();
      registerProtein_(protein_name);
      // score of group leader might technically not be transferable (due to
      // protein length etc.), but we still transfer it to allow filtering of
      // proteins by score without disrupting the groups:
      prot_id_->getHits().back().setScore(score);
    }
    else if (tag == "peptide")
    {
      // If a peptide is degenerate it will show in multiple groups, but have different statistics (e.g. 'nsp_adjusted_probability')
      // We thus treat each instance as a separate peptide
      // todo/improvement: link them by a group in PeptideIdentification?!
      pep_hit_ = new PeptideHit;
      pep_hit_->setSequence(AASequence::fromString(String(attributeAsString_(attributes, "peptide_sequence"))));
      pep_hit_->setScore(attributeAsDouble_(attributes, "nsp_adjusted_probability"));

      Int charge;
      if (optionalAttributeAsInt_(charge, attributes, "charge"))
      {
        pep_hit_->setCharge(charge);
      }
      else
      {
        OPENMS_LOG_WARN << "Required attribute 'charge' missing\n";
      }

      // add accessions of all indistinguishable proteins the peptide belongs to
      ProteinIdentification::ProteinGroup& indist = prot_id_->getIndistinguishableProteins().back();
      for (StringList::const_iterator accession = indist.accessions.begin(); accession != indist.accessions.end(); ++accession)
      {
        PeptideEvidence pe;
        pe.setProteinAccession(*accession);
        pep_hit_->addPeptideEvidence(pe);
      }
      pep_hit_->setMetaValue("is_unique", String(attributeAsString_(attributes, "is_nondegenerate_evidence")) == "Y" ? 1 : 0);
      pep_hit_->setMetaValue("is_contributing", String(attributeAsString_(attributes, "is_contributing_evidence")) == "Y" ? 1 : 0);
    }
    else if (tag == "mod_aminoacid_mass")
    {
      // relates to the last seen peptide (we hope)
      Size position = attributeAsInt_(attributes, "position");
      double mass = attributeAsDouble_(attributes, "mass");
      AASequence temp_aa_sequence(pep_hit_->getSequence());

      String temp_description = "";
      String origin = temp_aa_sequence[position - 1].getOneLetterCode();
      matchModification_(mass, origin, temp_description);
      if (!temp_description.empty()) // only if a mod was found
      {
        // e.g. Carboxymethyl (C)
        vector<String> mod_split;
        temp_description.split(' ', mod_split);
        if (mod_split.size() == 2)
        {
          if (mod_split[1] == "(C-term)" || ModificationsDB::getInstance()->getModification(temp_description)->getTermSpecificity() == ResidueModification::C_TERM)
          {
            temp_aa_sequence.setCTerminalModification(mod_split[0]);
          }
          else
          {
            if (mod_split[1] == "(N-term)" || ModificationsDB::getInstance()->getModification(temp_description)->getTermSpecificity() == ResidueModification::N_TERM)
            {
              temp_aa_sequence.setNTerminalModification(mod_split[0]);
            }
            else
            {
              // search this mod, if not directly use a general one
              temp_aa_sequence.setModification(position - 1, mod_split[0]);
            }
          }
        }
        else
        {
          error(LOAD, String("Cannot parse modification '") + temp_description + "@" + position + "'");
        }
      }
      else
      {
        error(LOAD, String("Cannot find modification '") + String(mass) + " " + String(origin) + "' @" + String(position));
      }

      pep_hit_->setSequence(temp_aa_sequence);
    }
  }

  void ProtXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    String tag = sm_.convert(qname);


    if (tag == "protein_group")
    {
      prot_id_->insertProteinGroup(protein_group_);
    }
    else if (tag == "peptide")
    {
      pep_id_->insertHit(*pep_hit_);
      delete pep_hit_;
    }
  }

  void ProtXMLFile::registerProtein_(const String& protein_name)
  {
    ProteinHit hit;
    hit.setAccession(protein_name);
    prot_id_->insertHit(hit);
    // add protein to groups
    protein_group_.accessions.push_back(protein_name);
    prot_id_->getIndistinguishableProteins().back().accessions.push_back(
      protein_name);
  }

  void ProtXMLFile::matchModification_(const double mass, const String& origin, String& modification_description)
  {
    double mod_mass = mass - ResidueDB::getInstance()->getResidue(origin)->getMonoWeight(Residue::Internal);
    vector<String> mods;
    ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, mod_mass, 0.001, origin);

    if (mods.size() == 1)
    {
      modification_description = mods[0];
    }
    else
    {
      if (!mods.empty())
      {
        String mod_str = mods[0];
        for (vector<String>::const_iterator mit = ++mods.begin(); mit != mods.end(); ++mit)
        {
          mod_str += ", " + *mit;
        }
        error(LOAD, "Modification '" + String(mass) + "' is not uniquely defined by the given data. Using '" + mods[0] +  "' to represent any of '" + mod_str + "'!");
        modification_description = mods[0];
      }
    }
  }

} // namespace OpenMS
