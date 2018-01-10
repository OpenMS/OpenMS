// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XQuestResultXMLHandler.h>
#include <xercesc/sax2/Attributes.hpp>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <boost/assign/list_of.hpp>
#include <iostream>
#include <utility>

#include <cassert>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
  namespace Internal
  {
    // Initialize static const members
    std::map< Size, String > XQuestResultXMLHandler::enzymes = boost::assign::map_list_of(0, "no_enzyme")
        (1, "trypsin") (2, "chymotrypsin") (3, "unknown_enzyme") (9, "unknown_enzyme")
        (10, "unknown_enzyme") (14, "unknown_enzyme") (15, "unknown_enzyme") (16, "unknown_enzyme") (17, "unknown_enzyme")
        (18, "unknown_enzyme") (20, "unknown_enzyme");

    std::map< String, UInt> XQuestResultXMLHandler::months = boost::assign::map_list_of("Jan", 1)
        ("Feb", 2) ("Mar", 3) ("Apr", 4) ("May", 5) ("Jun", 6) ("Jul", 7) ("Aug", 8) ("Sep", 9) ("Oct", 10) ("Nov", 11)("Dec", 12);

    const String XQuestResultXMLHandler::decoy_string = "decoy_";

    XQuestResultXMLHandler::XQuestResultXMLHandler(const String &filename,
                                                   std::vector< std::vector< PeptideIdentification > > & csms,
                                                   std::vector< ProteinIdentification > & prot_ids,
                                                   Size min_n_ions_per_spectrum,
                                                   bool load_to_peptideHit) :
      XMLHandler(filename, "1.0"),
      csms_(csms),
      prot_ids_(prot_ids),
      n_hits_(0),
      min_score_(0),
      max_score_(0),
      min_n_ions_per_spectrum_(min_n_ions_per_spectrum),
      load_to_peptideHit_(load_to_peptideHit)
    {
      // Initialize the one and only protein identification
      this->prot_ids_.clear();
      ProteinIdentification prot_id;
      prot_id.setSearchEngine("OpenXQuest");
      prot_id.setSearchEngineVersion(VersionInfo::getVersion());
      prot_id.setMetaValue("SpectrumIdentificationProtocol", DataValue("MS:1002494")); // cross-linking search = MS:1002494
      this->prot_ids_.push_back(prot_id);

      // Fetch the enzymes database
      this->enzymes_db_ = ProteaseDB::getInstance();

      // Produce some warnings that are associated with the reading of xQuest result files
      LOG_WARN << "WARNING: Fixed modifications are not available in the xQuest input file and will thus be not present in the loaded data!\n" << std::endl;
    }

    XQuestResultXMLHandler::~XQuestResultXMLHandler()
    {

    }

    void XQuestResultXMLHandler::extractDateTime_(const String & xquest_datetime_string, DateTime & date_time)
    {
      StringList xquest_datetime_string_split;
      StringUtils::split(xquest_datetime_string,' ', xquest_datetime_string_split);
      if (this->is_openproxl_)
      {
        // Example: 2017-03-17 23:04:50
        date_time.setDate(xquest_datetime_string_split[0]);
        date_time.setTime(xquest_datetime_string_split[1]);
      }
      else
      {
        // Example: Fri Dec 18 12:28:42 2015
        UInt day = xquest_datetime_string_split[2].toInt();
        UInt year = xquest_datetime_string_split[4].toInt();
        UInt month = XQuestResultXMLHandler::months[xquest_datetime_string_split[1]];
        date_time.setDate(month, day, year);
        date_time.setTime(xquest_datetime_string_split[3]);
      }
    }

    // Extracts the position of the Cross-Link for intralinks and crosslinks
    void XQuestResultXMLHandler::getLinkPosition_(const xercesc::Attributes & attributes, std::pair<SignedSize, SignedSize> & pair)
    {
      String xlink_position = this->attributeAsString_(attributes, "xlinkposition");
      StringList xlink_position_split;
      StringUtils::split(xlink_position, "," ,xlink_position_split);
      //assert(xlink_position_split.size() == 2);

      pair.first = xlink_position_split[0].toInt();
      pair.second = xlink_position_split.size() == 2 ? xlink_position_split[1].toInt() : 0;
    }

    void XQuestResultXMLHandler::setPeptideEvidence_(const String & prot_string, PeptideHit & pep_hit)
    {
      StringList prot_list;
      StringUtils::split(prot_string, ",", prot_list);
      vector< PeptideEvidence > evidences;
      evidences.reserve(prot_list.size());

      for (StringList::const_iterator prot_list_it = prot_list.begin();
           prot_list_it != prot_list.end(); ++prot_list_it)
      {
        PeptideEvidence pep_ev;
        String accession = *prot_list_it;

        if (this->accessions_.find(accession) == this->accessions_.end())
        {
          this->accessions_.insert(accession);

          ProteinHit prot_hit;
          prot_hit.setAccession(accession);
          prot_hit.setMetaValue("target_decoy", accession.hasSubstring("decoy") ? "decoy" : "target");

          this->prot_ids_[0].getHits().push_back(prot_hit);
        }

        pep_ev.setProteinAccession(accession);
        pep_ev.setStart(PeptideEvidence::UNKNOWN_POSITION); // These information are not available in the xQuest result file
        pep_ev.setEnd(PeptideEvidence::UNKNOWN_POSITION);
        pep_ev.setAABefore(PeptideEvidence::UNKNOWN_AA);
        pep_ev.setAAAfter(PeptideEvidence::UNKNOWN_AA);

        evidences.push_back(pep_ev);
      }
      pep_hit.setPeptideEvidences(evidences);
    }

    // Assign all attributes in the peptide_id_attributes map to the MetaInfoInterface object
    void XQuestResultXMLHandler::addMetaValues_(MetaInfoInterface & meta_info_interface)
    {
      for (std::map<String, DataValue>::const_iterator it = this->peptide_id_meta_values_.begin();
           it != this->peptide_id_meta_values_.end(); ++it)
      {
        std::pair<String, DataValue> item = *it;
        meta_info_interface.setMetaValue(item.first, item.second);
      }
    }

    void XQuestResultXMLHandler::setMetaValue_(const String &  key , const DataValue &  datavalue , PeptideIdentification & pep_id, PeptideHit & alpha)
    {
      pep_id.setMetaValue(key, datavalue);
      if (this->load_to_peptideHit_)
      {
        alpha.setMetaValue(key, datavalue);
      }
    }
    void XQuestResultXMLHandler::setMetaValue_(const String &  key , const DataValue &  datavalue , PeptideIdentification & pep_id, PeptideHit & alpha, PeptideHit & beta)
    {
      this->setMetaValue_(key, datavalue, pep_id, alpha);
      if (this->load_to_peptideHit_)
      {
        beta.setMetaValue(key, datavalue);
      }
    }

    double XQuestResultXMLHandler::getMinScore() const
    {
      return this->min_score_;
    }

    double XQuestResultXMLHandler::getMaxScore() const
    {
      return this->max_score_;
    }

    UInt XQuestResultXMLHandler::getNumberOfHits() const
    {
      return this->n_hits_;
    }

    void XQuestResultXMLHandler::endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname)
    {
      String tag = XMLString::transcode(qname);
      if (tag == "spectrum_search")
      {
        // Push back spectrum search vector
        Size current_spectrum_size = this->current_spectrum_search_.size();
        if (current_spectrum_size >= this->min_n_ions_per_spectrum_)
        {
          /* Currently the correct rank order within the xQuest file is assumed
            vector< PeptideIdentification > newvec(current_spectrum_size);
            for (vector< PeptideIdentification>::const_iterator it = this->current_spectrum_search.begin();
                it != this->current_spectrum_search.end(); ++it)
            {
              int index = (int) it->getMetaValue("xl_rank") - 1;

              if (newvec[index].metaValueExists("xl_rank"))
              {
                LOG_ERROR << "ERROR: At least two hits with the same rank within the same spectrum" << endl;
                throw std::exception();
              }
              newvec[index] = *it;
            }
            this->csms_.push_back(newvec);
            */
          this->csms_.push_back(this->current_spectrum_search_);
        }
        this->current_spectrum_search_.clear();
      }
      else if (tag == "xquest_results")
      {
        ProteinIdentification::SearchParameters search_params(this->prot_ids_[0].getSearchParameters());
        //search_params.charges = ListUtils::concatenate(this->charges_, ",");
        //search_params.setMetaValue("precursor:min_charge", this->min_precursor_charge_);
        //search_params.setMetaValue("precursor:max_charge", this->max_precursor_charge_);
        DoubleList monolink_masses_list(this->monolinks_masses_.size());
        for (std::set<double>::const_iterator monolink_masses_it = this->monolinks_masses_.begin();
             monolink_masses_it != this->monolinks_masses_.end(); ++monolink_masses_it)
        {
          monolink_masses_list.push_back(*monolink_masses_it);
        }
        search_params.setMetaValue("cross_link:mass_monolink", DataValue(monolink_masses_list));

        this->prot_ids_[0].setSearchParameters(search_params);
      }
    }
    void XQuestResultXMLHandler::startElement(const XMLCh * const, const XMLCh * const, const XMLCh * const qname, const Attributes &attributes)
    {
      String tag = XMLString::transcode(qname);
      // Extract meta information from the xquest_results tag
      if (tag == "xquest_results")
      {
        // Decide whether this Block is original xQuest or OpenProXL
        String xquest_version = this->attributeAsString_(attributes, "xquest_version");
        this->is_openproxl_ = xquest_version.hasSubstring("XL");

        // Date and Time of Search
        DateTime date_time;
        this->extractDateTime_(this->attributeAsString_(attributes, "date"), date_time);
        this->prot_ids_[0].setDateTime(date_time);

        // Set the search parameters
        ProteinIdentification::SearchParameters search_params;

        // General
        if (this->is_openproxl_) // Enzyme via name
        {
          search_params.digestion_enzyme = dynamic_cast<const DigestionEnzymeProtein&>(*this->enzymes_db_->getEnzyme(this->attributeAsString_(attributes, "enzyme_name")));
        }
        else // Enzyme via enzyme number in xQuest
        {
          search_params.digestion_enzyme = dynamic_cast<const DigestionEnzymeProtein&>(*this->enzymes_db_->getEnzyme(XQuestResultXMLHandler::enzymes[this->attributeAsInt_(attributes, "enzyme_num")]));
        }

        search_params.missed_cleavages = this->attributeAsInt_(attributes, "missed_cleavages");
        search_params.db = this->attributeAsString_(attributes, "database");
        search_params.precursor_mass_tolerance = this->attributeAsDouble_(attributes, "ms1tolerance");
        String tolerancemeasure = this->attributeAsString_(attributes, this->is_openproxl_ ? "tolerancemeasure_ms1" : "tolerancemeasure");
        search_params.precursor_mass_tolerance_ppm = tolerancemeasure == "ppm";
        search_params.fragment_mass_tolerance = this->attributeAsDouble_(attributes, "ms2tolerance");
        String tolerancemeasure_ms2 = this->attributeAsString_(attributes, "tolerancemeasure_ms2");
        search_params.fragment_mass_tolerance_ppm = tolerancemeasure_ms2 != "Da";

        // Modifications
        vector< String > variable_mod_list;
        vector< String > variable_mod_split;
        StringUtils::split(this->attributeAsString_(attributes, "variable_mod"), ",", variable_mod_split);

        // Oxidation of M
        if (variable_mod_split[0] == "M")
        {
          variable_mod_list.push_back("Oxidation (M)");
        }
        search_params.variable_modifications = variable_mod_list;

        // Meta Values
        search_params.setMetaValue("input_decoys", DataValue(this->attributeAsString_(attributes, "database_dc")));
        search_params.setMetaValue("decoy_prefix", DataValue(1));
        search_params.setMetaValue("decoy_string", DataValue(XQuestResultXMLHandler::decoy_string));
        search_params.setMetaValue("fragment:mass_tolerance_xlinks", DataValue(this->attributeAsDouble_(attributes, "xlink_ms2tolerance")));

        this->prot_ids_[0].setSearchParameters(search_params);
      }
      else if (tag == "spectrum_search")
      {

        // Update retention time of light
        StringList rt_split;
        StringUtils::split(this->attributeAsString_(attributes, "rtsecscans"), ":", rt_split);
        this->rt_light_ = rt_split[0].toDouble();

        // Update min and max precursor charge
        UInt charge_precursor = this->attributeAsInt_(attributes, "charge_precursor");
        if (charge_precursor < this->min_precursor_charge_)
        {
          this->min_precursor_charge_ = charge_precursor;
        }
        if (charge_precursor > this->max_precursor_charge_)
        {
          this->max_precursor_charge_ = charge_precursor;
        }
      }
      else if (tag == "search_hit")
      {
        // Keep track of the charge if this hit
        UInt charge = this->attributeAsInt_(attributes, "charge");
        this->charges_.insert(charge);

        this->n_hits_++;
        PeptideIdentification peptide_identification;

        // Set Attributes of Peptide Identification
        peptide_identification.setMZ(this->attributeAsDouble_(attributes, "mz"));
        peptide_identification.setRT(this->rt_light_);
        peptide_identification.setScoreType("OpenXQuest:combined score"); // Needed, since hard-coded in MzIdentMLHandler

        PeptideHit peptide_hit_alpha;
        PeptideHit peptide_hit_beta;
        vector<PeptideHit> peptide_hits;
        // XL Type, determined by "type"
        String xlink_type_string = this->attributeAsString_(attributes, "type");
        String prot1_string = this->attributeAsString_(attributes, "prot1");

        // Decide if decoy for alpha
        DataValue target_decoy = DataValue(prot1_string.hasSubstring("decoy") ? "decoy" : "target");
        peptide_identification.setMetaValue("target_decoy", target_decoy);
        peptide_hit_alpha.setMetaValue("target_decoy", target_decoy);

        // Set xl_chain meta value for alpha
        peptide_hit_alpha.setMetaValue("xl_chain", "MS:1002509");

        // Attributes of peptide_hit_alpha
        double score = this->attributeAsDouble_(attributes, "score");
        DataValue xlinkermass = DataValue(this->attributeAsDouble_(attributes, "xlinkermass"));

        // Set minscore and maxscore encountered
        if (score < this->min_score_)
        {
          this->min_score_ = score;
        }
        if (score > this->max_score_)
        {
          this->max_score_ = score;
        }
        peptide_hit_alpha.setScore(score);

        String seq1 = String(this->attributeAsString_(attributes, "seq1"));
        peptide_hit_alpha.setSequence(AASequence::fromString(seq1.substitute("X", "M(Oxidation)")));
        peptide_hit_alpha.setCharge(charge);

        // Get common attributes of Peptide Identification
        this->peptide_id_meta_values_["OpenXQuest:id"] = DataValue(this->attributeAsString_(attributes, "id"));
        this->peptide_id_meta_values_["OpenXQuest:xlinkermass"] = xlinkermass;
        this->peptide_id_meta_values_["OpenXQuest:wTIC"] = DataValue(this->attributeAsDouble_(attributes, "wTIC"));
        this->peptide_id_meta_values_["OpenXQuest:percTIC"] = DataValue(this->attributeAsDouble_(attributes, "TIC"));
        this->peptide_id_meta_values_["xl_rank"] = DataValue(this->attributeAsInt_(attributes, "search_hit_rank"));
        this->peptide_id_meta_values_["OpenXQuest:intsum"] = DataValue(this->attributeAsDouble_(attributes, "intsum")/100);
        this->peptide_id_meta_values_["OpenXQuest:match-odds"] = DataValue(this->attributeAsDouble_(attributes, "match_odds"));
        this->peptide_id_meta_values_["OpenXQuest:score"] = DataValue(score);
        this->peptide_id_meta_values_["OpenXQuest:error_rel"] = DataValue(this->attributeAsDouble_(attributes, "error_rel"));
        this->peptide_id_meta_values_["OpenXQuest:structure"] = DataValue(this->attributeAsString_(attributes, "structure"));

        assert(this->peptide_id_meta_values_["OpenXQuest:id"] != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_["OpenXQuest:xlinkermass"] != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_["OpenXQuest:wTIC"]  != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_["OpenXQuest:percTIC"] != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_["OpenXQuest:intsum"] != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_["OpenXQuest:match-odds"] != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_["xl_rank"] != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_["OpenXQuest:score"] != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_["OpenXQuest:error_rel"] != DataValue::EMPTY);
        assert(this->peptide_id_meta_values_["OpenXQuest:structure"] != DataValue::EMPTY);

        this->addMetaValues_(peptide_identification);

        // Store common attributes in Peptide Identification
        // If requested, also write to the peptide_hit_alpha
        if (this->load_to_peptideHit_)
        {
          this->addMetaValues_(peptide_hit_alpha);
        }

        // Store specific stuff for peptide hit alpha
        peptide_hit_alpha.setMetaValue("OpenXQuest:num_of_matched_ions",
                                       DataValue(this->attributeAsInt_(attributes, "num_of_matched_ions_alpha")));
        peptide_hit_alpha.setMetaValue("OpenXQuest:prot", DataValue(prot1_string));
        peptide_hit_alpha.setMetaValue("xl_mass", xlinkermass);

        // Set peptide Evidences for Alpha (need one for each accession in the prot1_string)
        this->setPeptideEvidence_(prot1_string, peptide_hit_alpha);

        // Switch on Cross-link type
        if (xlink_type_string == "xlink")
        {
          // Set the cross Link Mass
          ProteinIdentification::SearchParameters search_params(this->prot_ids_[0].getSearchParameters());
          search_params.setMetaValue("cross_link:mass", DataValue(this->attributeAsDouble_(attributes, "xlinkermass")));
          this->prot_ids_[0].setSearchParameters(search_params);

          peptide_hit_beta.setScore(score);

          String seq2 = String(this->attributeAsString_(attributes, "seq2"));
          peptide_hit_beta.setSequence(AASequence::fromString(seq2.substitute("X", "M(Oxidation)")));
          peptide_hit_beta.setCharge(charge);

          // If requested, also write to the peptide_hit_beta
          if (this->load_to_peptideHit_)
          {
            this->addMetaValues_(peptide_hit_beta);
          }
          // Set xl_type
          this->setMetaValue_("xl_type", DataValue("cross-link"), peptide_identification, peptide_hit_alpha, peptide_hit_beta);

          // Set xl positions, depends on xl_type
          std::pair<SignedSize, SignedSize> positions;
          this->getLinkPosition_(attributes, positions);
          peptide_hit_alpha.setMetaValue("xl_pos", DataValue(positions.first - 1));
          peptide_hit_beta.setMetaValue("xl_pos", DataValue(positions.second - 1));

          // Protein
          String prot2_string = this->attributeAsString_(attributes, "prot2");

          // Decide if decoy for beta
          if (prot2_string.hasSubstring("decoy"))
          {
            peptide_identification.setMetaValue("target_decoy", DataValue("decoy"));
            peptide_hit_beta.setMetaValue("target_decoy", DataValue("decoy"));
          }
          else
          {
            peptide_hit_beta.setMetaValue("target_decoy", DataValue("target"));
          }

          //  Set xl_chain meta value for beta
          peptide_hit_beta.setMetaValue("xl_chain", "MS:1002510");

          // Set peptide_hit specific stuff
          peptide_hit_beta.setMetaValue("OpenXQuest:num_of_matched_ions",
                                        DataValue(this->attributeAsInt_(attributes, "num_of_matched_ions_beta")));
          peptide_hit_beta.setMetaValue("OpenXQuest:prot", DataValue(prot2_string));
          peptide_hit_beta.setMetaValue("xl_mass", xlinkermass);

          // Set Peptide Evidences for Beta
          this->setPeptideEvidence_(prot2_string, peptide_hit_beta);

          // Determine if protein is intra/inter protein, check all protein ID combinations
          StringList prot1_list;
          prot1_string.split(",", prot1_list);
          StringList prot2_list;
          prot2_string.split( ",", prot2_list);

          for (StringList::const_iterator it1 = prot1_list.begin(); it1 != prot1_list.end(); ++it1)
          {
            for (StringList::const_iterator it2 = prot2_list.begin(); it2 != prot2_list.end(); ++it2)
            {
              String s1 = *it1;
              String s2 = *it2;
              s1.substitute("reverse_", "");
              s2.substitute("reverse_", "");
              s1.substitute(XQuestResultXMLHandler::decoy_string, "");
              s2.substitute(XQuestResultXMLHandler::decoy_string, "");

              this->setMetaValue_((s1.compare(s2) == 0) ? "OpenXQuest:is_intraprotein" : "OpenXQuest:is_interprotein",
                                 DataValue(), peptide_identification, peptide_hit_alpha, peptide_hit_beta);
            }
          }
        }
        else if (xlink_type_string == "intralink")
        {
          // xl type
          this->setMetaValue_("xl_type", DataValue("loop-link"), peptide_identification, peptide_hit_alpha);

          // Set xl positions, depends on xl_type
          std::pair<SignedSize, SignedSize> positions;
          this->getLinkPosition_(attributes, positions);
          peptide_hit_alpha.setMetaValue("xl_pos", DataValue(positions.first - 1));
          peptide_hit_alpha.setMetaValue("xl_pos2", DataValue(positions.second - 1));
        }
        else if (xlink_type_string == "monolink")
        {
          // Set the monolink mass
          this->monolinks_masses_.insert(this->attributeAsDouble_(attributes, "xlinkermass"));

          // xl_type
          this->setMetaValue_("xl_type", DataValue("mono-link"),peptide_identification,peptide_hit_alpha);

          std::pair< SignedSize, SignedSize > xlink_pos;
          this->getLinkPosition_(attributes, xlink_pos);
          peptide_hit_alpha.setMetaValue("xl_pos", DataValue(xlink_pos.first - 1));
        }
        else
        {
          LOG_ERROR << "ERROR: Unsupported Cross-Link type: " << xlink_type_string << endl;
          throw std::exception();
        }

        // Finalize this record
        peptide_hits.push_back(peptide_hit_alpha);

        if (peptide_hit_beta.metaValueExists("xl_pos"))
        {
          peptide_hits.push_back(peptide_hit_beta);
        }

        peptide_identification.setHits(peptide_hits);
        this->peptide_id_meta_values_.clear();
        this->current_spectrum_search_.push_back(peptide_identification);
      }
    }
  }   // namespace Internal
} // namespace OpenMS
