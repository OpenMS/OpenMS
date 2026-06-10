// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace xercesc;
using namespace std;

namespace OpenMS
{

  XTandemXMLFile::XTandemXMLFile() :
    XMLHandler("", "1.1"),
    XMLFile()
  {
    // see X! Tandem parameters "protein, quick pyrolidone" and "protein, quick acetyl":
    default_nterm_mods_.setModifications("", "Gln->pyro-Glu (N-term Q),Glu->pyro-Glu (N-term E),Acetyl (N-term)");
  }

  XTandemXMLFile::~XTandemXMLFile() = default;

  void XTandemXMLFile::load(const std::string& filename, ProteinIdentification& protein_identification, PeptideIdentificationList& peptide_ids, ModificationDefinitionsSet& mod_def_set)
  {
    // File name for error message in XMLHandler
    file_ = filename;
    mod_def_set_ = mod_def_set;

    // reset everything, in case "load" is called multiple times:
    is_protein_note_ = is_spectrum_note_ = skip_protein_acc_update_ = false;
    peptide_hits_.clear();
    protein_hits_.clear();
    current_protein_ = tag_ = previous_seq_ = "";
    current_charge_ = 0;
    current_id_ = current_start_ = current_stop_ = 0;
    spectrum_ids_.clear();

    enforceEncoding_("ISO-8859-1");
    parse_(filename, this);

    DateTime now = DateTime::now();
    std::string date_string = now.getDate();
    std::string identifier("XTandem_" + date_string);
    //vector<std::string> accessions;

    // convert mapping id -> peptide_hits into peptide hits list
    peptide_ids.clear();
    for (map<UInt, vector<PeptideHit> >::iterator it = peptide_hits_.begin(); it != peptide_hits_.end(); ++it)
    {
      PeptideIdentification id;
      id.setScoreType("XTandem");
      id.setHigherScoreBetter(true);
      id.setIdentifier(identifier);
      id.setSpectrumReference( spectrum_ids_[it->first]);

      id.getHits().swap(it->second);
      id.sort();
      peptide_ids.push_back(id);
    }

    protein_identification.getHits().swap(protein_hits_);

    // E-values
    protein_identification.setHigherScoreBetter(false);
    protein_identification.sort();
    protein_identification.setScoreType("XTandem");
    protein_identification.setSearchEngine("XTandem");

    // TODO version of XTandem ???? is not available from performance param section of outputfile (to be parsed)
    // TODO Date of search, dito
    protein_identification.setDateTime(now);
    protein_identification.setIdentifier(identifier);

    // TODO search parameters are also available

    // mods may be changed, copy them back:
    mod_def_set = mod_def_set_;
  }

  void XTandemXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
  {
    tag_ =std::string(sm_.convert(qname));

    if (tag_ == "domain")
    {
      std::string id_string = attributeAsString_(attributes, "id");
      UInt id = StringUtils::toInt32(StringUtils::prefix(id_string, '.'));
      current_id_ = id;

      PeptideEvidence pe;
      // get amino acid before
      std::string pre = attributeAsString_(attributes, "pre");
      if (!pre.empty())
      {
        pe.setAABefore(pre[pre.size()-1]);
      }
      // get amino acid after
      std::string post = attributeAsString_(attributes, "post");
      if (!post.empty())
      {
        pe.setAAAfter(post[0]);
      }

      current_start_ = attributeAsInt_(attributes, "start");
      pe.setStart(current_start_ - 1);
      current_stop_ = attributeAsInt_(attributes, "end");
      pe.setEnd(current_stop_ - 1);
        
      pe.setProteinAccession(current_protein_);
 
      std::string seq = attributeAsString_(attributes, "seq");
      // is this the same peptide as before, just in a different protein (scores will be the same)?
      if ((!peptide_hits_.contains(id)) || (seq != previous_seq_))
      {
        PeptideHit hit;
        // can't parse sequences permissively because that would skip characters
        // that X! Tandem includes when calculating e.g. modification positions,
        // potentially leading to errors when assigning mods to residues:
        hit.setSequence(AASequence::fromString(seq, false));
        hit.setCharge(current_charge_);
        
        // get scores etc.:
        hit.setMetaValue("nextscore", attributeAsDouble_(attributes, "nextscore"));
        hit.setMetaValue("delta", attributeAsDouble_(attributes, "delta"));
        hit.setMetaValue("mass", attributeAsDouble_(attributes, "mh")); // note the different names
        hit.setMetaValue("E-Value", attributeAsDouble_(attributes, "expect")); // note the different names
        double hyperscore = attributeAsDouble_(attributes, "hyperscore");
        hit.setMetaValue("hyperscore", hyperscore);
        hit.setScore(hyperscore);

        // try to get a, b, c, x, y, z score (optional)
        std::string ions = "abcxyz";
        std::string ion_score = " _score";
        std::string ion_count = " _ions";
        for (auto it = ions.begin(); it != ions.end(); ++it)
        {
          ion_score[0] = *it;
          double score;
          if (optionalAttributeAsDouble_(score, attributes, ion_score.c_str()))
          {
            hit.setMetaValue(ion_score, score);
          }
          ion_count[0] = *it;
          UInt count;
          if (optionalAttributeAsUInt_(count, attributes, ion_count.c_str()))
          {
            hit.setMetaValue(ion_count, count);
          }
        }

        peptide_hits_[id].push_back(hit);
        previous_seq_ = seq;
      }

      peptide_hits_[id].back().addPeptideEvidence(pe);
      return;
    }

    if (tag_ == "aa")
    {
      if (group_type_stack_.empty())
      {
        error(LOAD,std::string("Found an 'aa' element outside of a 'group'! Please check your input file."));
      }
      auto& current_group_type_ = group_type_stack_.top();
      // TODO support "aa" entries in the parameter groups (e.g. to read user-specified amino acids)
      //  currently we just ignore them
      if (current_group_type_ == GroupType::MODEL)
      {
        // e.g. <aa type="S" at="2" modified="42.0106" />
        std::string aa = attributeAsString_(attributes, "type");
        Int mod_pos = attributeAsInt_(attributes, "at");
        double mass_shift = attributeAsDouble_(attributes, "modified");

        AASequence aa_seq = peptide_hits_[current_id_].back().getSequence();
        mod_pos -= current_start_; // X! Tandem uses position in the protein

        const ResidueModification* res_mod = nullptr;

        // first, try to find matching mod in the defined ones:
        multimap<double, ModificationDefinition> matches;
        if (mod_pos <= 0) // N-terminal mod?
        {
          mod_def_set_.findMatches(matches, mass_shift, aa,
                                   ResidueModification::N_TERM);
        }
        else if (mod_pos >= Int(aa_seq.size() - 1)) // C-terminal mod?
        {
        mod_def_set_.findMatches(matches, mass_shift, aa,
                                    ResidueModification::C_TERM);
        }
        if (matches.empty())
        {
          mod_def_set_.findMatches(matches, mass_shift, aa,
                                   ResidueModification::ANYWHERE);
        }
        if (matches.empty() && (mod_pos <= 0)) // try X! Tandem's default mods
        {
          default_nterm_mods_.findMatches(matches, mass_shift, aa,
                                          ResidueModification::N_TERM);
          if (!matches.empty())
          {
            // add the match to the mod. set for output (search parameters):
            ModificationDefinition mod_def(matches.begin()->second.getModificationName());
            mod_def.setFixedModification(false);
            mod_def_set_.addModification(mod_def);
          }
        }
        // matches are sorted by mass error - first one is best match:
        if (!matches.empty())
        {
          res_mod = &(matches.begin()->second.getModification());
        }
        else // no match? strange, let's try all possible modifications
        {
          const ModificationsDB* mod_db = ModificationsDB::getInstance();
          // try to find a mod that fits
          if (mod_pos <= 0) // can (!) be an N-terminal mod
          {
            res_mod = mod_db->getBestModificationByDiffMonoMass(mass_shift, 0.01, aa, ResidueModification::N_TERM);
          }
          else if (mod_pos >= static_cast<int>(aa_seq.size() - 1))
            {
              res_mod = mod_db->getBestModificationByDiffMonoMass(mass_shift, 0.01, aa, ResidueModification::C_TERM);
            }
          if (res_mod == nullptr) // if no terminal mod, try normal one
          {
            res_mod = mod_db->getBestModificationByDiffMonoMass(mass_shift, 0.01, aa, ResidueModification::ANYWHERE);
          }
        }

        if (res_mod == nullptr)
        {
          error(LOAD,std::string("No modification found which fits residue '") + aa + "' with mass '" + StringUtils::toStr(mass_shift) + "'!");
        }
        else
        {
          // @TODO: avoid unnecessary conversion of mods from/to string below
          if (res_mod->getTermSpecificity() == ResidueModification::N_TERM)
          {
            aa_seq.setNTerminalModification(res_mod->getFullId());
          }
          else if (res_mod->getTermSpecificity() == ResidueModification::C_TERM)
          {
            aa_seq.setCTerminalModification(res_mod->getFullId());
          }
          else
          {
            aa_seq.setModification(mod_pos, res_mod->getFullId());
          }
          peptide_hits_[current_id_].back().setSequence(aa_seq);
        }
      }
    }

    if (tag_ == "group")
    {
      std::string type = attributeAsString_(attributes, "type");
      if (type == "model")
      {
        group_type_stack_.push(GroupType::MODEL);
        Int index = attributes.getIndex(sm_.convert("z").c_str());
        if (index >= 0)
        {
          current_charge_ = StringUtils::toInt32(std::string(sm_.convert(attributes.getValue(index))));
        }
        previous_seq_ = "";
      }
      else if (type == "parameter")
      {
        group_type_stack_.push(GroupType::PARAMETERS);
      }
      else
      {
        group_type_stack_.push(GroupType::SUPPORT);
      }
    }

    if (tag_ == "note")
    {
      std::string label;
      optionalAttributeAsString_(label, attributes, "label");

      if (label == "description") // in '<"protein" ...>'
      {
        is_protein_note_ = true;
      }
      else if (label == "Description") // in '<group type="support" label="fragment ion mass spectrum">'
      {
        is_spectrum_note_ = true;
      }
    }

    if (tag_ == "protein")
    {
      UInt uid = attributeAsInt_(attributes, "uid");
      if (!protein_uids_.contains(uid)) // new protein
      {
        ProteinHit hit;
        // accession may be overwritten based on '<note label="description">', but set it for now:
        current_protein_ = attributeAsString_(attributes, "label");
        hit.setAccession(current_protein_);

        double score(0);
        if (optionalAttributeAsDouble_(score, attributes, "expect"))
        {
          hit.setScore(score);
        }
        protein_hits_.push_back(hit);
        protein_uids_.insert(uid);
        skip_protein_acc_update_ = false;
      }
      else
      {
        current_protein_ = attributeAsString_(attributes, "label"); //save it for upcoming peptides (domains)
        skip_protein_acc_update_ = true;
      }
      return;
    }

  }

  void XTandemXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    tag_ =std::string(sm_.convert(qname));
    if (tag_ == "group")
    {
      group_type_stack_.pop();
    }
  }

  void XTandemXMLFile::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
  {
    if (tag_ == "note")
    {
      if (is_protein_note_)
      {
        current_protein_ =StringUtils::trimmed(std::string(sm_.convert(chars)));
        if (!skip_protein_acc_update_)
        {
          protein_hits_.back().setAccession(current_protein_);
        }
      }
      else if (is_spectrum_note_)
      {
        spectrum_ids_[current_id_] =StringUtils::trimmed(std::string(sm_.convert(chars)));
      }
      is_protein_note_ = false;
      is_spectrum_note_ = false;
    }
  }

} // namespace OpenMS
