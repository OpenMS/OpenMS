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
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XQuestResultXMLHandler.h>
#include <xercesc/sax2/Attributes.hpp>
#include <OpenMS/METADATA/XQuestResultMeta.h>
#include <iostream>
#include <utility>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
  namespace Internal
  {

    XQuestResultXMLHandler::XQuestResultXMLHandler(const String &filename,
                                                   vector< XQuestResultMeta> & metas,
                                                   std::vector< std::vector< CrossLinkSpectrumMatch > > & csms,
                                                   int & n_hits_,
                                                   std::vector< int > * cum_hits) :
      XMLHandler(filename, "1.0"),
      metas_(metas),
      csms_(csms),
      n_hits_(n_hits_),
      cum_hits_(cum_hits)
    {
    }
    XQuestResultXMLHandler::~XQuestResultXMLHandler()
    {

    }

    void XQuestResultXMLHandler::endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname)
    {
      String tag = XMLString::transcode(qname);

      if (tag == "spectrum_search")
      {
          // Push back spectrum search vector and ensure that the hits are sorted by their rank within the vector
          vector< CrossLinkSpectrumMatch > newvec(this->current_spectrum_search.size());
          for(vector< CrossLinkSpectrumMatch>::const_iterator it = this->current_spectrum_search.begin();
              it != this->current_spectrum_search.end(); ++it)
          {
            if (newvec[it->rank - 1].rank != 0)
            {
              LOG_ERROR << "ERROR: At least two hits with the same rank within the same spectrum" << endl;
              throw std::exception();
            }
            newvec[it->rank - 1] = *it;
          }
          this->csms_.push_back(newvec);
          this->current_spectrum_search.clear();
          // Might be NULL if we do not want to calculate the cum_hits in the first place
          if(this->cum_hits_ != NULL)
          {
            this->cum_hits_->push_back(this->n_hits_);\
          }
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
          CrossLinkSpectrumMatch csm;
          ProteinProteinCrossLink cross_link;

          // XL Type, determined by "type"
          ProteinProteinCrossLink::ProteinProteinCrossLinkType xlink_type;
          String xlink_type_string = this->attributeAsString_(attributes, "type");

          if (xlink_type_string == "xlink")
          {

              xlink_type = ProteinProteinCrossLink::CROSS;
          }
          else if (xlink_type_string == "monolink")
          {

              xlink_type = ProteinProteinCrossLink::MONO;
          }
          else if (xlink_type_string == "intralink")
          {

              xlink_type = ProteinProteinCrossLink::LOOP;
          }
          else
          {
            LOG_ERROR << "ERROR: Unsupported Cross-Link type: " << xlink_type_string << endl;
            throw std::exception();

          }


          // Switch on the type of XLINK
          if (xlink_type == ProteinProteinCrossLink::CROSS)
          {

            AASequence seq1 = AASequence::fromString(this->attributeAsString_(attributes, "seq1"));
            AASequence seq2 = AASequence::fromString(this->attributeAsString_(attributes, "seq2"));

            // Ensure that alpha is not shorter than beta
            if (seq1.size() > seq2.size())
            {
              cross_link.alpha = seq1;
              cross_link.beta = seq2;
            }
            else
            {
              cross_link.alpha = seq2;
              cross_link.beta = seq1;
            }
          }
          else
          {
            cross_link.alpha = AASequence::fromString(this->attributeAsString_(attributes, "seq1"));
            assert(cross_link.beta.empty());
          }

          if (xlink_type == ProteinProteinCrossLink::MONO)
          {
              // Important: (cross_link_position.second == -1)
              cross_link.cross_link_position = std::make_pair(this->attributeAsInt_(attributes, "xlinkposition"), -1);
          }
          else
          {
              String xlink_position = this->attributeAsString_(attributes, "xlinkposition");
              StringList xlink_position_split;
              StringUtils::split(xlink_position, "," ,xlink_position_split);
              if (xlink_position_split.size() != 2)
              {
                String message = "ERROR: Nonsense specification of cross-link position: " + xlink_position;
                LOG_ERROR << message << endl;
                throw std::exception();
              }
              cross_link.cross_link_position = std::make_pair(xlink_position_split[0].toInt(),
                                                              xlink_position_split[1].toInt());
          }

          cross_link.cross_linker_mass = this->attributeAsDouble_(attributes, "xlinkermass");
          cross_link.cross_linker_name = "NA";   // Cross-Linker name cannot be recovered from XQuest Result XML

          csm.cross_link = cross_link;
          csm.percTIC = this->attributeAsDouble_(attributes, "TIC");
          csm.wTIC = this->attributeAsDouble_(attributes, "wTIC");
          csm.int_sum = this->attributeAsDouble_(attributes, "intsum") / 100;
          csm.match_odds = this->attributeAsDouble_(attributes, "match_odds");
          csm.xcorrx_max = this->attributeAsDouble_(attributes, "xcorrx");
          csm.xcorrc_max = this->attributeAsDouble_(attributes, "xcorrb"); // ?????????
          csm.rank = this->attributeAsInt_(attributes, "search_hit_rank");
          csm.score = this->attributeAsDouble_(attributes, "score");
          csm.error_rel = this->attributeAsDouble_(attributes, "error_rel");
          csm.structure = this->attributeAsString_(attributes, "structure");
          csm.num_of_matched_ions_alpha = this->attributeAsInt_(attributes, "num_of_matched_ions_alpha");

          // Missing number for beta is indicated with '-'
          try
          {
            csm.num_of_matched_ions_beta = this->attributeAsInt_(attributes, "num_of_matched_ions_beta");
          }
          catch (...)
          {
            csm.num_of_matched_ions_beta  = -1;
          }
          this->current_spectrum_search.push_back(csm);
      }
    }
    void XQuestResultXMLHandler::characters(const XMLCh * const chars, const XMLSize_t)
    {
    }


  }   // namespace Internal
} // namespace OpenMS










