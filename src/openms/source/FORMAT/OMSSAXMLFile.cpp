// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSSAXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

using namespace std;

namespace OpenMS
{

  OMSSAXMLFile::OMSSAXMLFile() :
    XMLHandler("", 1.1),
    XMLFile(),
    peptide_identifications_(nullptr)
  {
    readMappingFile_();
  }

  OMSSAXMLFile::~OMSSAXMLFile() = default;

  void OMSSAXMLFile::load(const String& filename, ProteinIdentification& protein_identification, vector<PeptideIdentification>& peptide_identifications, bool load_proteins, bool load_empty_hits)
  {
    // clear input (in case load() is called more than once)
    protein_identification = ProteinIdentification();
    peptide_identifications.clear();

    // filename for error messages in XMLHandler
    file_ = filename;

    load_proteins_ = load_proteins;
    load_empty_hits_ = load_empty_hits;
    peptide_identifications_ = &peptide_identifications;

    // fill the internal datastructures
    parse_(filename, this);

    DateTime now = DateTime::now();
    String identifier("OMSSA_" + now.get());

    // post-processing
    set<String> accessions;
    for (PeptideIdentification& pep : peptide_identifications)
    {
      pep.setScoreType("OMSSA");
      pep.setHigherScoreBetter(false);
      pep.setIdentifier(identifier);
      pep.assignRanks();

      if (load_proteins)
      {
        for (const PeptideHit& pit : pep.getHits())
        {
          set<String> hit_accessions = pit.extractProteinAccessionsSet();
          accessions.insert(hit_accessions.begin(), hit_accessions.end());
        }
      }
    }

    if (load_proteins)
    {
      for (set<String>::const_iterator it = accessions.begin(); it != accessions.end(); ++it)
      {
        ProteinHit hit;
        hit.setAccession(*it);
        protein_identification.insertHit(hit);
      }

      // E-values
      protein_identification.setHigherScoreBetter(false);
      protein_identification.setScoreType("OMSSA");
      protein_identification.setIdentifier(identifier);
    }

    // version of OMSSA is not available
    // Date of the search is not available -> set it to now()
    protein_identification.setDateTime(now);
    protein_identification.setIdentifier(identifier);

    // search parameters are also not available
  }

  void OMSSAXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& /*attributes*/)
  {
    tag_ = String(sm_.convert(qname)).trim();
  }

  void OMSSAXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    tag_ = String(sm_.convert(qname)).trim();

    // protein hits (MSPepHits) are handled in ::characters(...)

    // end of peptide hit
    if (tag_ == "MSHits")
    {
      actual_peptide_hit_.setPeptideEvidences(actual_peptide_evidences_);
      actual_peptide_evidence_ = PeptideEvidence();
      actual_peptide_evidences_.clear();
      actual_peptide_id_.insertHit(actual_peptide_hit_);
      actual_peptide_hit_ = PeptideHit();
    }
    // end of peptide id
    else if (tag_ == "MSHitSet")
    {
      if (!actual_peptide_id_.getHits().empty()  || load_empty_hits_)
      {
        peptide_identifications_->push_back(actual_peptide_id_);
      }
      actual_peptide_id_ = PeptideIdentification();
    }
    else if (tag_ == "MSModHit")
    {
      /*
        Modifications:
        <MSModHit>
          <MSModHit_site>1</MSModHit_site>
          <MSModHit_modtype>
            <MSMod>3</MSMod>
          </MSModHit_modtype>
        </MSModHit>
      */
      if (mods_map_.find(actual_mod_type_.toInt()) != mods_map_.end() && !mods_map_[actual_mod_type_.toInt()].empty())
      {
        if (mods_map_[actual_mod_type_.toInt()].size() > 1)
        {
          warning(LOAD, String("Cannot determine exact type of modification of position ") + actual_mod_site_ + " in sequence " + actual_peptide_hit_.getSequence().toString() + " using modification " + actual_mod_type_ + " - using first possibility!");
        }
        AASequence pep = actual_peptide_hit_.getSequence();
        auto mod = *(mods_map_[actual_mod_type_.toInt()].begin());
        if (mod->getTermSpecificity() == ResidueModification::N_TERM)
        {
          pep.setNTerminalModification(mod->getFullId());
        }
        else if (mod->getTermSpecificity() == ResidueModification::C_TERM)
        {
          pep.setCTerminalModification(mod->getFullId());
        }
        else
        {
          pep.setModification(actual_mod_site_, mod->getFullId());
        }
        actual_peptide_hit_.setSequence(pep);
      }
      else
      {
        warning(LOAD, String("Cannot find PSI-MOD mapping for mod - ignoring '") + actual_mod_type_ + "'");
      }
    }

    tag_ = "";
  }

  void OMSSAXMLFile::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
  {
    if (tag_.empty()) return;

    String value = ((String)sm_.convert(chars)).trim();
    // MSPepHit section
    // <MSPepHit_start>0</MSPepHit_start>
    // <MSPepHit_stop>8</MSPepHit_stop>
    // <MSPepHit_accession>6599</MSPepHit_accession>
    // <MSPepHit_defline>CRHU2 carbonate dehydratase (EC 4.2.1.1) II [validated] - human</MSPepHit_defline>
    // <MSPepHit_protlength>260</MSPepHit_protlength>
    // <MSPepHit_oid>6599</MSPepHit_oid>

    if (tag_ == "MSPepHit_start")
    {
      tag_ = "";
      return;
    }
    else if (tag_ == "MSPepHit_stop")
    {
      tag_ = "";
      return;
    }
    else if (tag_ == "MSPepHit_accession")
    {
      if (load_proteins_)
      {
        actual_peptide_evidence_.setProteinAccession(value);
      }
      tag_ = "";
      return;
    }
    else if (tag_ == "MSPepHit_defline")
    {
      // TODO add defline to ProteinHit?
      tag_ = "";
      return;
    }
    else if (tag_ == "MSPepHit_protlength")
    {
      tag_ = "";
      return;
    }
    else if (tag_ == "MSPepHit_oid")
    {
      tag_ = "";
      // end of MSPepHit so we add the evidence
      actual_peptide_evidences_.push_back(actual_peptide_evidence_);
      return;
    }
    // MSHits section
    // <MSHits_evalue>0.00336753988893542</MSHits_evalue>
    // <MSHits_pvalue>1.30819399070598e-08</MSHits_pvalue>
    // <MSHits_charge>1</MSHits_charge>
    // <MSHits_pepstring>MSHHWGYGK</MSHits_pepstring>
    // <MSHits_mass>1101492</MSHits_mass>
    // <MSHits_pepstart></MSHits_pepstart>
    // <MSHits_pepstop>H</MSHits_pepstop>
    // <MSHits_theomass>1101484</MSHits_theomass>
    else if (tag_ == "MSHits_evalue")
    {
      actual_peptide_hit_.setScore(value.toDouble());
      tag_ = "";
      return;
    }
    else if (tag_ == "MSHits_charge")
    {
      actual_peptide_hit_.setCharge(value.toInt());
      tag_ = "";
      return;
    }
    else if (tag_ == "MSHits_pvalue")
    {
      // TODO extra field?
      //actual_peptide_hit_.setScore(value.toDouble());
      tag_ = "";
      return;
    }
    else if (tag_ == "MSHits_pepstring")
    {
      AASequence seq;
      try
      {
        seq = AASequence::fromString(value.trim());
      }
      catch (Exception::ParseError& /* e */)
      {
        actual_peptide_hit_.setSequence(AASequence());
        tag_ = "";
        return;
      }

      if (mod_def_set_.getNumberOfFixedModifications() != 0)
      {
        set<String> fixed_mod_names = mod_def_set_.getFixedModificationNames();
        for (set<String>::const_iterator it = fixed_mod_names.begin(); it != fixed_mod_names.end(); ++it)
        {
          String origin = ModificationsDB::getInstance()->getModification(*it)->getOrigin();
          UInt position(0);
          for (const Residue& ait : seq)
          {
            if (ait.getOneLetterCode() == origin)
            {
              seq.setModification(position, *it);
            }
            ++position;
          }
        }
      }
      actual_peptide_hit_.setSequence(seq);
      tag_ = "";
      return;
    }
    else if (tag_ == "MSHits_mass")
    {
      tag_ = "";
      return;
    }
    else if (tag_ == "MSHits_pepstart")
    {
      if (!value.empty() && !actual_peptide_evidences_.empty())
      {
        actual_peptide_evidences_[0].setAABefore(value[0]);
      }
      tag_ = "";
      return;
    }
    else if (tag_ == "MSHits_pepstop")
    {
      if (!value.empty() && !actual_peptide_evidences_.empty())
      {
        actual_peptide_evidences_[0].setAAAfter(value[0]);
      }
      tag_ = "";
      return;
    }
    else if (tag_ == "MSHits_theomass")
    {
      tag_ = "";
      return;
    }

    // modifications
    ///<MSHits_mods>
    // <MSModHit>
    //  <MSModHit_site>4</MSModHit_site>
    //  <MSModHit_modtype>
    //   <MSMod>2</MSMod>
    //  </MSModHit_modtype>
    // </MSModHit>
    ///</MSHits_mods>


    if (tag_ == "MSHits_mods")
    {
      actual_mod_site_ = 0;
      actual_mod_type_ = "";
    }
    else if (tag_ == "MSModHit_site")
    {
      actual_mod_site_ = value.trim().toInt();
    }
    else if (tag_ == "MSMod")
    {
      actual_mod_type_ = value.trim();
    }
    // m/z value and rt
    else if (tag_ == "MSHitSet_ids_E")
    {
      // value might be  ( OMSSA 2.1.8): 359.213256835938_3000.13720000002_controllerType=0 controllerNumber=1 scan=4655
      //                 (<OMSSA 2.1.8): 359.213256835938_3000.13720000002
      if (!value.trim().empty())
      {
        if (value.has('_'))
        {
          StringList sp = ListUtils::create<String>(value, '_');
          try
          {
            actual_peptide_id_.setMZ(sp[0].toDouble());
            actual_peptide_id_.setRT(sp[1].toDouble());
          }
          catch (...)
          {
            // if exception happens to occur here, s.th. went wrong, e.g. the value does not contain numbers
          }
        }
      }
    }
  }

  void OMSSAXMLFile::readMappingFile_()
  {
    String file = File::find("CHEMISTRY/OMSSA_modification_mapping");
    TextFile infile(file);

    for (TextFile::ConstIterator it = infile.begin(); it != infile.end(); ++it)
    {
      vector<String> split;
      it->split(',', split);

      if (!it->empty() && (*it)[0] != '#')
      {
        Int omssa_mod_num = split[0].trim().toInt();
        if (split.size() < 2)
        {
          fatalError(LOAD, String("Invalid mapping file line: '") + *it + "'");
        }
        vector<const ResidueModification*> mods;
        for (Size i = 2; i != split.size(); ++i)
        {
          String tmp(split[i].trim());
          if (!tmp.empty())
          {
            const ResidueModification* mod = ModificationsDB::getInstance()->getModification(tmp);
            mods.push_back(mod);
            mods_to_num_[mod->getFullId()] = omssa_mod_num;
          }
        }
        mods_map_[omssa_mod_num] = mods;
      }
    }
  }

  void OMSSAXMLFile::setModificationDefinitionsSet(const ModificationDefinitionsSet& mod_set)
  {
    mod_def_set_ = mod_set;
    UInt omssa_mod_num(119);
    set<String> mod_names = mod_set.getVariableModificationNames();
    for (set<String>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
    {
      if (!(mods_to_num_.find(*it) != mods_to_num_.end()))
      {
        mods_map_[omssa_mod_num].push_back(ModificationsDB::getInstance()->getModification(*it));
        mods_to_num_[*it] = omssa_mod_num;
        ++omssa_mod_num;
      }
    }
  }

} // namespace OpenMS
