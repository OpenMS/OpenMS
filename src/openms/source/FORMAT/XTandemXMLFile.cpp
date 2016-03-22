// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace xercesc;
using namespace std;

namespace OpenMS
{

  XTandemXMLFile::XTandemXMLFile() :
    XMLHandler("", 1.1),
    XMLFile(),
    actual_start_(0),
    actual_stop_(0)
  {

  }

  XTandemXMLFile::~XTandemXMLFile()
  {
  }

  void XTandemXMLFile::setModificationDefinitionsSet(const ModificationDefinitionsSet& rhs)
  {
    mod_def_set_ = rhs;
  }

  void XTandemXMLFile::load(const String& filename, ProteinIdentification& protein_identification, vector<PeptideIdentification>& peptide_ids)
  {
    //File name for error message in XMLHandler
    file_ = filename;

    enforceEncoding_("ISO-8859-1");
    parse_(filename, this);

    DateTime now = DateTime::now();
    String date_string = now.getDate();
    String identifier("XTandem_" + date_string);
    //vector<String> accessions;

    // convert id -> peptide_hits into peptide hits list
    //vector<PeptideIdentification> peptide_identifications;
    PeptideIdentification().metaRegistry().registerName("spectrum_id", "the id of the spectrum counting from 1");
    for (map<UInt, vector<PeptideHit> >::const_iterator it = peptide_hits_.begin(); it != peptide_hits_.end(); ++it)
    {
      // reduce the hits with the same sequence to one PeptideHit
      map<String, vector<PeptideHit> > seq_to_hits;
      for (vector<PeptideHit>::const_iterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
      {
        seq_to_hits[it1->getSequence().toString()].push_back(*it1);
      }

      PeptideIdentification id;
      // if (descriptions_.find(it->first) != descriptions_.end())
      // {
      // id.setMetaValue("Description", descriptions_[it->first]);
      // }
      for (map<String, vector<PeptideHit> >::const_iterator it1 = seq_to_hits.begin(); it1 != seq_to_hits.end(); ++it1)
      {
        const vector<PeptideHit>& peptide_hits = it->second;
        if (!peptide_hits.empty())
        {
          // store all peptide hits that identify the same sequence in a single peptide hit
          PeptideHit hit = peptide_hits[0];
          vector<PeptideEvidence> peptide_evidences;
          for (vector<PeptideHit>::const_iterator it2 = peptide_hits.begin(); it2 != peptide_hits.end(); ++it2)
          {
            const vector<PeptideEvidence> evidences = it2->getPeptideEvidences();
            for (vector<PeptideEvidence>::const_iterator e_it = evidences.begin(); e_it != evidences.end(); ++e_it)
            {
              // only rewrite accession and keep AABefore/AAAfter, start, stop information from peptide evidence
              PeptideEvidence pe = *e_it;
              String new_acc = protein_hits_[pe.getProteinAccession()].getAccession();
              pe.setProteinAccession(new_acc);
              peptide_evidences.push_back(pe);
            }
          }
          hit.setPeptideEvidences(peptide_evidences);
          id.insertHit(hit);
        }
      }

      id.setScoreType("XTandem");
      id.setHigherScoreBetter(true);
      id.setIdentifier(identifier);
      id.assignRanks();
      id.setMetaValue("spectrum_id", it->first);

      peptide_ids.push_back(id);
    }

    //sort(accessions.begin(), accessions.end());
    //vector<String>::const_iterator end_unique = unique(accessions.begin(), accessions.end());

    for (Map<String, ProteinHit>::const_iterator pit = protein_hits_.begin(); pit != protein_hits_.end(); ++pit)
    {
      protein_identification.insertHit(pit->second);
    }


    // E-values
    protein_identification.setHigherScoreBetter(false);
    protein_identification.assignRanks();
    protein_identification.setScoreType("XTandem");
    protein_identification.setSearchEngine("XTandem");

    // TODO version of XTandem ???? is not available from performance param section of outputfile (to be parsed)
    // TODO Date of search, dito
    protein_identification.setDateTime(now);
    protein_identification.setIdentifier(identifier);


    // TODO search parameters are also available
  }

  void XTandemXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
  {
    tag_ = String(sm_.convert(qname));

    if (tag_ == "domain")
    {
      PeptideEvidence pe;
      PeptideHit hit;
      hit.metaRegistry().registerName("E-Value", "E-Value of hit");
      hit.metaRegistry().registerName("nextscore", "next_score of hit");
      // get nextscore
      double nextscore(String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("nextscore"))))).toDouble());
      hit.setMetaValue("nextscore", nextscore);
      // get hyperscore
      double hyperscore(String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("hyperscore"))))).toDouble());
      hit.setScore(hyperscore);

      // get mass
      double mass(String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("mh"))))).toDouble());
      hit.setMetaValue("mass", mass);

      // get delta
      double delta(String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("delta"))))).toDouble());
      hit.setMetaValue("delta", delta);

      // try to get a, b, c, x, y, z score (optional)
      String att_str;
      if (optionalAttributeAsString_(att_str, attributes, "a_score")) hit.setMetaValue("a_score", att_str.toDouble());
      if (optionalAttributeAsString_(att_str, attributes, "a_ions")) hit.setMetaValue("a_ions", att_str.toInt());
      
      if (optionalAttributeAsString_(att_str, attributes, "b_score")) hit.setMetaValue("b_score", att_str.toDouble());
      if (optionalAttributeAsString_(att_str, attributes, "b_ions")) hit.setMetaValue("b_ions", att_str.toInt());

      if (optionalAttributeAsString_(att_str, attributes, "c_score")) hit.setMetaValue("c_score", att_str.toDouble());
      if (optionalAttributeAsString_(att_str, attributes, "c_ions")) hit.setMetaValue("c_ions", att_str.toInt());

      if (optionalAttributeAsString_(att_str, attributes, "x_score")) hit.setMetaValue("x_score", att_str.toDouble());
      if (optionalAttributeAsString_(att_str, attributes, "x_ions")) hit.setMetaValue("x_ions", att_str.toInt());

      if (optionalAttributeAsString_(att_str, attributes, "y_score")) hit.setMetaValue("y_score", att_str.toDouble());
      if (optionalAttributeAsString_(att_str, attributes, "y_ions")) hit.setMetaValue("y_ions", att_str.toInt());

      if (optionalAttributeAsString_(att_str, attributes, "z_score")) hit.setMetaValue("z_score", att_str.toDouble());
      if (optionalAttributeAsString_(att_str, attributes, "z_ions")) hit.setMetaValue("z_ions", att_str.toInt());

      // get sequence of peptide
      String seq = attributeAsString_(attributes, "seq");
      hit.setSequence(AASequence::fromString(seq));

      // get amino acid before
      String pre = attributeAsString_(attributes, "pre");
      if (!pre.empty())
      {
        pe.setAABefore(pre[pre.size()-1]);
      }

      // get amino acid after
      String post = attributeAsString_(attributes, "post");
      if (!post.empty())
      {
        pe.setAAAfter(post[0]);
      }

      // get expectation value
      String expect = attributeAsString_(attributes, "expect");
      hit.setMetaValue("E-Value", expect.toDouble());

      // get precursor m/z
      //double mh(String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("mh"))))).toDouble());
      //hit.setMetaValue("MZ", mh); // not needed, set by the XTandem Adapter itself

      // spectrum id
      String id_string = attributeAsString_(attributes, "id");
      UInt id(id_string.prefix('.').toInt());
      // hit.setMetaValue("RT_index", id);
      actual_id_ = id;

      String tmp;
      optionalAttributeAsString_(tmp, attributes, "start");
      actual_start_ = tmp.toInt();
      pe.setStart(actual_start_ - 1);
      tmp = "";
      optionalAttributeAsString_(tmp, attributes, "end");
      actual_stop_ = tmp.toInt();
      pe.setEnd(actual_stop_ - 1);

      // add the actual protein accession
      pe.setProteinAccession(actual_protein_id_);
      hit.setCharge(actual_charge_);

      hit.addPeptideEvidence(pe);
      peptide_hits_[actual_id_].push_back(hit);
      return;
    }

    if (tag_ == "aa")
    {
      // e.g. <aa type="S" at="2" modified="42.0106" />
      String type, at, modified;
      optionalAttributeAsString_(type, attributes, "type");
      optionalAttributeAsString_(at, attributes, "at");
      optionalAttributeAsString_(modified, attributes, "modified");

      AASequence aa_seq = peptide_hits_[actual_id_].back().getSequence();
      UInt mod_pos = (UInt)at.toInt() - actual_start_;

      // search mod
      vector<String> possible_mods, possible_mass_mods;

      // try to find a mod in the given mods that fits

      if (mod_pos == 0) // can (!) be a N-terminal mod
      {
        ModificationsDB::getInstance()->getTerminalModificationsByDiffMonoMass(possible_mass_mods, modified.toDouble(), 0.01, ResidueModification::N_TERM);
      }
      else if (mod_pos == aa_seq.size())
      {
        ModificationsDB::getInstance()->getTerminalModificationsByDiffMonoMass(possible_mass_mods, modified.toDouble(), 0.01, ResidueModification::C_TERM);
      }

      // if not found a terminal mod, try normal one
      if (possible_mass_mods.empty())
      {
        ModificationsDB::getInstance()->getModificationsByDiffMonoMass(possible_mass_mods, type, modified.toDouble(), 0.01);
      }

      // cerr << "Possible mods of type='" << type << "', weight='" << modified.toDouble() << "', mod_pos='" << mod_pos << "'" << "\n";
      // for (vector<String>::const_iterator it = possible_mass_mods.begin(); it != possible_mass_mods.end(); ++it)
      // {
      // cerr << *it << " " << ModificationsDB::getInstance()->getModification(*it).getTermSpecificity() << "\n";
      // }

      set<String> mod_names = mod_def_set_.getModificationNames();

      // throw out any of the modifications that are not contained in the def set (throws out also s.th. like "Carbamidomethyl (N-term)"
      for (vector<String>::const_iterator it = possible_mass_mods.begin(); it != possible_mass_mods.end(); ++it)
      {
        if (mod_names.find(*it) != mod_names.end())
        {
          possible_mods.push_back(*it);
        }
      }

      // cerr << "Possible mods (#=" << possible_mods.size() << "): " << "\n";
      // for (vector<String>::const_iterator it = possible_mods.begin(); it != possible_mods.end(); ++it)
      // {
      // cerr << *it << "\n";
      // }

      // maybe we missed the real modification, even if it's not terminal
      if (possible_mods.empty() && mod_pos == 0)
      {
        vector<String> new_possible_mass_mods;
        ModificationsDB::getInstance()->getModificationsByDiffMonoMass(new_possible_mass_mods, type, modified.toDouble(), 0.01);
        // now try to find this in the definitions which are set
        for (vector<String>::const_iterator it = new_possible_mass_mods.begin(); it != new_possible_mass_mods.end(); ++it)
        {
          if (mod_names.find(*it) != mod_names.end())
          {
            possible_mods.push_back(*it);
          }
        }
      }


      // use all possible mass mods, because the modification was not predefined
      if (possible_mods.empty())
      {
        possible_mods = possible_mass_mods;
      }

      if (possible_mods.empty())
      {
        error(LOAD, String("No modification found which fits residue '") + type + "' with mass '" + modified + "'!");
      }
      else
      {
        if (possible_mods.size() > 1)
        {
          // if available use a specific one, except if the modification is terminal
          set<String> specific_ones;
          for (vector<String>::const_iterator it = possible_mods.begin(); it != possible_mods.end(); ++it)
          {
            String origin  = ModificationsDB::getInstance()->getModification(*it).getOrigin();
            ResidueModification::Term_Specificity term_spec = ModificationsDB::getInstance()->getModification(*it).getTermSpecificity();
            if (origin == type && !(term_spec == ResidueModification::N_TERM || term_spec == ResidueModification::C_TERM))
            {
              specific_ones.insert(*it);
            }
          }

          if (specific_ones.size() == 1)
          {
            possible_mods.clear();
            possible_mods.push_back(*specific_ones.begin());
          }
          else
          {
            if (specific_ones.empty())
            {
              // maybe there are terminal modifications but none of them has been selected
              // search unspecific but terminal residues
              vector<String> new_possible_mods;
              for (vector<String>::const_iterator it = possible_mods.begin(); it != possible_mods.end(); ++it)
              {
                String origin  = ModificationsDB::getInstance()->getModification(*it).getOrigin();
                ResidueModification::Term_Specificity term_spec = ModificationsDB::getInstance()->getModification(*it).getTermSpecificity();
                //cerr << "Testing: " << *it << ", origin='" << origin << "', term_spec='" << term_spec << "'" << "\n";
                if ((origin == "N-term" || origin == "C-term") && (term_spec == ResidueModification::N_TERM || term_spec == ResidueModification::C_TERM))
                {
                  //cerr << "Adding1: '" << *it << "', origin='" << origin << "' term_spec='" << term_spec << "'" << "\n";
                  new_possible_mods.push_back(*it);
                }
              }

              if (new_possible_mods.empty())
              {
                // if we haven't found a generic terminal modification, we search for a specific terminal mods which fits
                for (vector<String>::const_iterator it = possible_mods.begin(); it != possible_mods.end(); ++it)
                {
                  String origin  = ModificationsDB::getInstance()->getModification(*it).getOrigin();
                  ResidueModification::Term_Specificity term_spec = ModificationsDB::getInstance()->getModification(*it).getTermSpecificity();
                  if (origin == type && (term_spec == ResidueModification::N_TERM || term_spec == ResidueModification::C_TERM))
                  {
                    //cerr << "Adding2: '" << *it << "', origin='" << origin << "' term_spec='" << term_spec << "'" << "\n";
                    new_possible_mods.push_back(*it);
                  }
                }
              }
              if (!new_possible_mods.empty())
              {
                possible_mods = new_possible_mods;
              }
            }
            else
            {
              // also haven't found a non-specific terminal modification
              //if (specific_ones.empty())
              //{
              // put the specific ones in front of the list
              vector<String> new_possible_mods;
              for (set<String>::const_iterator it = specific_ones.begin(); it != specific_ones.end(); ++it)
              {
                new_possible_mods.push_back(*it);
              }
              for (vector<String>::const_iterator it = possible_mods.begin(); it != possible_mods.end(); ++it)
              {
                if (specific_ones.find(*it) == specific_ones.end())
                {
                  new_possible_mods.push_back(*it);
                }
              }
              possible_mods = new_possible_mods;
              // }
              // else
              // {
              // for (set<String>::const_iterator it = specific_ones.begin(); it != specific_ones.end(); ++it)
              // {
              // possible_mods.push_back(*it);
              // }
            }
          }
          if (possible_mods.size() > 1)
          {
            String error_string = String("More than one modification found which fits residue '") + type + "' with mass '" + modified + "': ";
            String possbile_mods;
            possbile_mods.concatenate(possible_mods.begin(), possible_mods.end(), ',');
            error_string += possbile_mods + ". Using first hit: '" + *possible_mods.begin() + "'.";
            error(LOAD, error_string);
          }
        }

        if (ModificationsDB::getInstance()->getModification(*possible_mods.begin()).getTermSpecificity() == ResidueModification::N_TERM && mod_pos == 0)
        {
          aa_seq.setNTerminalModification(*possible_mods.begin());
        }
        else
        {
          aa_seq.setModification(mod_pos, *possible_mods.begin());
        }
      }

      peptide_hits_[actual_id_].back().setSequence(AASequence(aa_seq));

      return;
    }

    if (tag_ == "group")
    {
      Int index = attributes.getIndex(sm_.convert("z"));
      if (index >= 0)
      {
        actual_charge_ = String(sm_.convert(attributes.getValue(index))).toInt();
      }
      return;
    }

    if (tag_ == "note")
    {
      String label;
      optionalAttributeAsString_(label, attributes, "label");

      if (label == "description")
      {
        is_description_ = true;
      }
    }

    if (tag_ == "protein")
    {
      //protein_open_ = true;
      ProteinHit hit;

      String uid;
      optionalAttributeAsString_(uid, attributes, "uid");
      actual_protein_id_ = uid;

      if (!protein_hits_.has(uid))
      {
        double score(0);
        optionalAttributeAsDouble_(score, attributes, "expect");
        hit.setScore(score);

        protein_hits_[uid] = hit;
      }
      return;
    }

  }

  void XTandemXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    tag_ = String(sm_.convert(qname));
    return;
  }

  void XTandemXMLFile::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
  {
    if (tag_ == "note" && is_description_)
    {
      is_description_ = false;
      protein_hits_[actual_protein_id_].setAccession(((String) sm_.convert(chars)).trim());
    }
  }

} // namespace OpenMS
