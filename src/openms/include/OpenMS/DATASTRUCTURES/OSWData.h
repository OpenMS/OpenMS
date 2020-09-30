// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>
#include <vector>

namespace OpenMS
{
    struct OSWTransition
    {
      public:
        OSWTransition(const String& annotation, const UInt32 id , const float product_mz, const char type, const bool is_decoy)
          : annotation_(annotation),
            id_(id),
            product_mz_(product_mz),
            type_(type),
            is_decoy_(is_decoy)
        {}

        OSWTransition(const OSWTransition& rhs) = default;
        OSWTransition& operator=(const OSWTransition& rhs) = default;
        OSWTransition(OSWTransition&& rhs) = default;
        OSWTransition& operator=(OSWTransition&& rhs) = default;
        ~OSWTransition() = default;

        const String& getAnnotation() const
        {
          return annotation_;
        }
        UInt32 getID() const
        {
          return id_;
        }
        float getProductMZ() const
        {
          return product_mz_;
        }
        char getType() const
        {
          return type_;
        }
        bool isDecoy() const
        {
          return is_decoy_;
        }

      private:
        String annotation_; ///< e.g. y5/-0.002
        UInt32 id_;         /// ID as used in OSWPeakGroup::transition_ids
        float product_mz_;
        char type_;         ///< b, y,
        bool is_decoy_;
    };

    struct OSWPeakGroup // Feature == group of transtions in certain RT range
    {
      OSWPeakGroup(const float rt_experimental, const float rt_left_width, const float rt_right_width, const float rt_delta, const std::vector<UInt32>& transition_ids, const float q_value = -1)
        : rt_experimental_(rt_experimental),
          rt_left_width_(rt_left_width),
          rt_right_width_(rt_right_width),
          rt_delta_(rt_delta),
          q_value_(q_value),
          transition_ids_(transition_ids)
      {
      }
      float rt_experimental_;   ///< rt apex of this feature in seconds (averaged across all transitions)
      float rt_left_width_;     ///< rt start in seconds
      float rt_right_width_;    ///< rt end in seconds
      float rt_delta_;          ///< rt offset from expected distance
      float q_value_ = -1;      ///< optional Q-value from pyProphet; equals -1 if not set;
      std::vector<UInt32> transition_ids_; /// many features will point to the same transition (but at different RT);
    };

    struct OSWPeptidePrecursor
    {

      OSWPeptidePrecursor(const String& seq, const short charge, const bool decoy, const float precursor_mz, const std::vector<OSWPeakGroup>& features)
        : seq_(seq_),
          charge_(charge),
          decoy_(decoy),
          precursor_mz_(precursor_mz),
          features_(features)
      {
      }

      String seq_;
      short charge_;
      bool decoy_;
      float precursor_mz_;
      std::vector<OSWPeakGroup> features_;
    };

    struct OSWProtein
    {
      OSWProtein(const String& accession, const std::vector<OSWPeptidePrecursor>& peptides)
        : accession_(accession),
          peptides_(peptides)
      {}

      String accession_;
      std::vector<OSWPeptidePrecursor> peptides_;
    };

    class OSWData
    {
      public:

        /// Adds a transition; do this before adding Proteins
        void addTransition(const OSWTransition& tr)
        {
          transitions_.emplace(tr.getID(), tr);
        }

        /// Adds a protein, which has all its subcomponents already populated
        /// The transitions references internally are checked to make sure
        /// they are valid. 
        /// @throws Exception::Precondition() if transition IDs are unknown
        void addProtein(OSWProtein&& prot)
        {
          // check if transitions are known
          for (const auto& pc : prot.peptides_)
          {
            for (const auto& f : pc.features_)
            {
              for (const auto& tr : f.transition_ids_)
              {
                if (transitions_.find(tr) == transitions_.end())
                {
                  throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Transition with ID " + String(tr) + " was referenced in Protein/Precursor/Feature but is not known!");
                }
              }
            }
          }
          proteins_.push_back(std::move(prot));
        }

        std::vector<OSWProtein>::const_iterator protBegin() const
        {
          return proteins_.cbegin();
        }

        std::vector<OSWProtein>::const_iterator protEnd() const
        {
          return proteins_.cend();
        }

        Size protCount() const
        {
          return proteins_.size();
        }

        Size transitionCount() const
        {
          return transitions_.size();
        }

        const OSWTransition& getTransition(const UInt32 id) const
        {
          return transitions_.at(id);
        }
        
        std::map<UInt32, OSWTransition>::const_iterator transitionsBegin() const
        {
          return transitions_.cbegin();
        }

        std::map<UInt32, OSWTransition>::const_iterator transitionsEnd() const
        {
          return transitions_.cend();
        }

        /// forget all data
        void clear()
        {
          transitions_.clear();
          proteins_.clear();
        }

      private:

      std::map<UInt32, OSWTransition> transitions_;
      std::vector<OSWProtein> proteins_;
    };
    

} // namespace OpenMS