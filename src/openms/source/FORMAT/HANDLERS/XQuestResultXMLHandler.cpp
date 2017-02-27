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
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XQuestResultXMLHandler.h>
#include <xercesc/sax2/Attributes.hpp>
#include <OpenMS/METADATA/XQuestResultMeta.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <iostream>
#include <utility>

using namespace std;
using namespace xercesc;

// TODO Switch to Debug mode in CMake and remove this undef
#undef NDEBUG
#include <assert.h>

namespace OpenMS
{
  namespace Internal
  {

    /*
     *  Helper functions
     */
    void removeSubstring(String &  large, const String & small)
    {
      std::string::size_type i = large.find(small);
      if (i != std::string::npos)
      {
        large.erase(i, small.length());
      }
    }

    XQuestResultXMLHandler::XQuestResultXMLHandler(const String &filename,
                                                   vector< XQuestResultMeta> & metas,
                                                   std::vector< std::vector< PeptideIdentification > > & csms,
                                                   int & n_hits_,
                                                   std::vector< int > * cum_hits,
                                                   size_t min_n_ions_per_spectrum,
                                                   bool load_to_peptideHit) :
      XMLHandler(filename, "1.0"),
      metas_(metas),
      csms_(csms),
      n_hits_(n_hits_),
      cum_hits_(cum_hits),
      min_n_ions_per_spectrum_(min_n_ions_per_spectrum),
      load_to_peptideHit_(load_to_peptideHit)
    {
    }
    XQuestResultXMLHandler::~XQuestResultXMLHandler()
    {

    }


    // Extracts the position of the Cross-Link for intralinks and crosslinks
    void XQuestResultXMLHandler::getLinkPosition_(const xercesc::Attributes & attributes, std::pair<SignedSize, SignedSize> & pair)
    {
      String xlink_position = this->attributeAsString_(attributes, "xlinkposition");
      StringList xlink_position_split;
      StringUtils::split(xlink_position, "," ,xlink_position_split);
      assert(xlink_position_split.size() == 2);

      pair.first = xlink_position_split[0].toInt();
      pair.second = xlink_position_split[1].toInt();
    }


    void XQuestResultXMLHandler::endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname)
    {
      String tag = XMLString::transcode(qname);

      if (tag == "spectrum_search")
      {
          // Push back spectrum search vector and ensure that the hits are sorted by their rank within the vector
          size_t current_spectrum_size = this->current_spectrum_search.size();
          if (current_spectrum_size >= this->min_n_ions_per_spectrum_)
          {
            vector< PeptideIdentification > newvec(current_spectrum_size);
            for(vector< PeptideIdentification>::const_iterator it = this->current_spectrum_search.begin();
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
            if(this->cum_hits_ != NULL)
            {
              this->cum_hits_->push_back(this->n_hits_);\
            }
        }
        this->current_spectrum_search.clear();
      }
      else if (tag == "xquest_results")
      {
          this->metas_.push_back(this->current_meta_);
          this->current_meta_.clearMetaInfo();
      }
    }
    void XQuestResultXMLHandler::startElement(const XMLCh * const, const XMLCh * const, const XMLCh * const qname, const Attributes &attributes)
    {
      String tag = XMLString::transcode(qname);

      // Extract meta information
      if (tag == "xquest_results")
      {
        // For now, put each of the key value pair as DataValue into the meta interface
        for(XMLSize_t i = 0; i < attributes.getLength(); ++i)
        {
            this->current_meta_.setMetaValue(XMLString::transcode(attributes.getQName(i)),
                                     DataValue(XMLString::transcode(attributes.getValue(i))));
        }
      }
      else if (tag == "spectrum_search")
      {
        // TODO Store information of spectrum search
      }
      else if (tag == "search_hit")
      {
          // New CrossLinkSpectrumMatchEntry
          this->n_hits_++;
          PeptideIdentification peptide_identification;
          PeptideHit peptide_hit_alpha;
          vector<PeptideHit> peptide_hits;
          // XL Type, determined by "type"
          String xlink_type_string = this->attributeAsString_(attributes, "type");
          //String prot1 = this->attributeAsString_(attributes, "prot1");


          /* TODO Determine Decoy and intra/inter cross links
          if (prot1.hasSubstring("decoy"))
          {
              peptide_hit_alpha.setMetaValue("OpenXQuest:is_decoy", DataValue());
          }
          // That is a bit hacky. Maybe use a regular expressions instead.
          removeSubstring(prot1, "reverse_");
          removeSubstring(prot1, "decoy_");
          */

          // Get Attributes of Peptide Identification
          DataValue xquest_xlinkermass = DataValue(this->attributeAsDouble_(attributes, "xlinkermass"));
          DataValue xquest_wtic = DataValue(this->attributeAsDouble_(attributes, "wTIC"));
          DataValue xquest_perctic = DataValue(this->attributeAsDouble_(attributes, "TIC"));
          DataValue xquest_intsum = DataValue(this->attributeAsDouble_(attributes, "intsum")/100);
          DataValue xquest_match_odds = DataValue(this->attributeAsDouble_(attributes, "match_odds"));
          DataValue xquest_xl_rank = DataValue(this->attributeAsInt_(attributes, "search_hit_rank"));
          DataValue xquest_score = DataValue(this->attributeAsDouble_(attributes, "score"));
          DataValue xquest_error_rel = DataValue(this->attributeAsDouble_(attributes, "error_rel"));
          DataValue structure = DataValue(this->attributeAsString_(attributes, "structure"));

          assert(xquest_xlinkermass != DataValue::EMPTY);
          assert(xquest_wtic  != DataValue::EMPTY);
          assert(xquest_perctic != DataValue::EMPTY);
          assert(xquest_intsum != DataValue::EMPTY);
          assert(xquest_match_odds != DataValue::EMPTY);
          assert(xquest_xl_rank != DataValue::EMPTY);
          assert(xquest_score != DataValue::EMPTY);
          assert(xquest_error_rel != DataValue::EMPTY);
          assert(structure != DataValue::EMPTY);


          // Store attributes in Peptide Identification
          peptide_identification.setMetaValue("xl_rank", xquest_xl_rank);
          peptide_identification.setMetaValue("OpenXQuest:xlinkermass", xquest_xlinkermass);
          peptide_identification.setMetaValue("OpenXQuest:wTIC", xquest_wtic);
          peptide_identification.setMetaValue("OpenXQuest:percTIC", xquest_perctic);
          peptide_identification.setMetaValue("OpenXQuest:intsum", xquest_intsum);
          peptide_identification.setMetaValue("OpenXQuest:match-odds", xquest_match_odds);
          peptide_identification.setMetaValue("OpenXQuest:score", xquest_score);
          peptide_identification.setMetaValue("OpenXQuest:error_rel", xquest_error_rel);
          peptide_identification.setMetaValue("OpenXQuest:structure", structure);

          // If requested, also write to the peptide_hit_alpha
          if (this->load_to_peptideHit_)
          {
              peptide_hit_alpha.setMetaValue("xl_rank", xquest_xl_rank);
              peptide_hit_alpha.setMetaValue("OpenXQuest:xlinkermass", xquest_xlinkermass);
              peptide_hit_alpha.setMetaValue("OpenXQuest:wTIC", xquest_wtic);
              peptide_hit_alpha.setMetaValue("OpenXQuest:percTIC", xquest_perctic);
              peptide_hit_alpha.setMetaValue("OpenXQuest:intsum", xquest_intsum);
              peptide_hit_alpha.setMetaValue("OpenXQuest:match-odds", xquest_match_odds);
              peptide_hit_alpha.setMetaValue("OpenXQuest:score", xquest_score);
              peptide_hit_alpha.setMetaValue("OpenXQuest:error_rel", xquest_error_rel);
              peptide_hit_alpha.setMetaValue("OpenXQuest:structure", structure);
          }

          // Figure out cross-link type
          if (xlink_type_string == "xlink")
          {
              PeptideHit peptide_hit_beta;

              // If requested, also write to the peptide_hit_beta
              if (this->load_to_peptideHit_)
              {
                  peptide_hit_beta.setMetaValue("xl_rank", xquest_xl_rank);
                  peptide_hit_beta.setMetaValue("OpenXQuest:xlinkermass", xquest_xlinkermass);
                  peptide_hit_beta.setMetaValue("OpenXQuest:wTIC", xquest_wtic);
                  peptide_hit_beta.setMetaValue("OpenXQuest:percTIC", xquest_perctic);
                  peptide_hit_beta.setMetaValue("OpenXQuest:intsum", xquest_intsum);
                  peptide_hit_beta.setMetaValue("OpenXQuest:match-odds", xquest_match_odds);
                  peptide_hit_beta.setMetaValue("OpenXQuest:score", xquest_score);
                  peptide_hit_beta.setMetaValue("OpenXQuest:error_rel", xquest_error_rel);
                  peptide_hit_beta.setMetaValue("OpenXQuest:structure", structure);
              }

              peptide_identification.setMetaValue("xl_type", "cross-link");
              std::pair<SignedSize, SignedSize> positions;
              this->getLinkPosition_(attributes, positions);
              peptide_hit_alpha.setMetaValue("xl_pos", DataValue(positions.first));
              peptide_hit_beta.setMetaValue("xl_pos", DataValue(positions.second));

              // Number of matched ions
              peptide_hit_alpha.setMetaValue("OpenXQuest:num_of_matched_ions",
                                             DataValue(this->attributeAsInt_(attributes, "num_of_matched_ions_alpha")));
              peptide_hit_beta.setMetaValue("OpenXQuest:num_of_matched_ions",
                                            DataValue(this->attributeAsInt_(attributes, "num_of_matched_ions_beta")));
              if (this->load_to_peptideHit_)
              {
                peptide_hit_alpha.setMetaValue("xl_type", "cross-link");
                peptide_hit_beta.setMetaValue("xl_type", "cross-link");
              }
              peptide_hits.push_back(peptide_hit_beta);
          }
          else if (xlink_type_string == "intralink")
          {
              peptide_identification.setMetaValue("xl_type", "loop-link");

              std::pair<SignedSize, SignedSize> positions;
              this->getLinkPosition_(attributes, positions);
              peptide_hit_alpha.setMetaValue("xl_pos", DataValue(positions.first));
              peptide_hit_alpha.setMetaValue("xl_pos2", DataValue(positions.second));

              if (this->load_to_peptideHit_)
              {
                peptide_hit_alpha.setMetaValue("xl_type", "loop-link");
              }

              // Number of matched ions
              peptide_hit_alpha.setMetaValue("OpenXQuest:num_of_matched_ions",
                                             DataValue(this->attributeAsInt_(attributes, "num_of_matched_ions_alpha")));

          }
          else if (xlink_type_string == "monolink")
          {
             peptide_identification.setMetaValue("xl_type", "mono-link");
             peptide_hit_alpha.setMetaValue("xl_pos", DataValue((SignedSize)this->attributeAsInt_(attributes, "xlinkposition")));

             if (this->load_to_peptideHit_)
             {
               peptide_hit_alpha.setMetaValue("xl_type", "mono-link");
             }
             // Number of matched ions
             peptide_hit_alpha.setMetaValue("OpenXQuest:num_of_matched_ions",
                                            DataValue(this->attributeAsInt_(attributes, "num_of_matched_ions_alpha")));

          }
          else
          {
            LOG_ERROR << "ERROR: Unsupported Cross-Link type: " << xlink_type_string << endl;
            throw std::exception();
          }

          peptide_hits.push_back(peptide_hit_alpha);
          peptide_identification.setHits(peptide_hits);
          this->current_spectrum_search.push_back(peptide_identification);
      }
    }
    void XQuestResultXMLHandler::characters(const XMLCh * const chars, const XMLSize_t)
    {
    }


  }   // namespace Internal
} // namespace OpenMS










