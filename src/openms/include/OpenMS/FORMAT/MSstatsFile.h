// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Lukas Heumos $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <boost/regex.hpp>

namespace OpenMS
{
  /**
    @brief File adapter for MzTab files
    @ingroup FileIO
  */
  
    class OPENMS_DLLAPI MSstatsFile
    {
    public:
        ///Default constructor
        MSstatsFile();
        ///Destructor
        ~MSstatsFile();

        /// store label free experiment (MSstats)
        void storeLFQ(const String& filename, 
                      const ConsensusMap &consensus_map,
                      const ExperimentalDesign& design,
                      const StringList& reannotate_filenames,
                      const bool is_isotope_label_type,
                      const String& bioreplicate,
                      const String& condition,
                      const String& retention_time_summarization_method);
        
        /// store isobaric experiment (MSstatsTMT)
        void storeISO(const String& filename, 
                      const ConsensusMap &consensus_map,
                      const ExperimentalDesign& design,
                      const StringList& reannotate_filenames,
                      const String& bioreplicate,
                      const String& condition,
                      const String& mixture,
                      const String& retention_time_summarization_method);
    
    private:
     
        const String na_string = "NA";
        
        /// The meta value of the peptide identification which is going to be used for the experimental design link
        const String meta_value_exp_design_key = "spectra_data";

        /*
        *  @brief: Internal function to check if MSstats_BioReplicate and MSstats_Condition in Experimental Design
        */
        static void checkConditionLFQ_(const ExperimentalDesign::SampleSection& sampleSection, const String& bioreplicate, const String& condition);

        /*
         *  @brief: Internal function to check if MSstats_BioReplicate, MSstats_Condition and MSstats_Mixture in Experimental Design
         */
        static void checkConditionISO_(const ExperimentalDesign::SampleSection sampleSection, const String& bioreplicate, const String& condition, const String& mixture);

        /*
         *  @brief MSstats treats runs differently than OpenMS. In MSstats, runs are an enumeration of (SpectraFilePath, Fraction)
         *  In OpenMS, a run is split into multiple fractions.
         *
         */
        static void assembleRunMap(
                std::map< std::pair< String, unsigned>, unsigned> &run_map,
                const ExperimentalDesign &design)
        {
          run_map.clear();
          const ExperimentalDesign::MSFileSection& msfile_section = design.getMSFileSection();
          unsigned run_counter = 1;

          for (ExperimentalDesign::MSFileSectionEntry const& r : msfile_section)
          {
            std::pair< String, unsigned> tpl = std::make_pair(File::basename(r.path), r.fraction);
            if (run_map.find(tpl) == run_map.end())
            {
              run_map[tpl] = run_counter++;
            }
          }
        }

        bool checkUnorderedContent_(const std::vector< String> &first, const std::vector< String > &second)
        {
          const std::set< String > lhs(first.begin(), first.end());
          const std::set< String > rhs(second.begin(), second.end());
          return lhs == rhs
                 && std::equal(lhs.begin(), lhs.end(), rhs.begin());
        }

        OpenMS::Peak2D::IntensityType sumIntensity(const std::set< OpenMS::Peak2D::IntensityType > &intensities)
        {
          OpenMS::Peak2D::IntensityType result = 0;
          for (const OpenMS::Peak2D::IntensityType &intensity : intensities)
          {
            result += intensity;
          }
          return result;
        }

        OpenMS::Peak2D::IntensityType meanIntensity(const std::set< OpenMS::Peak2D::IntensityType > &intensities)
        {
          return sumIntensity(intensities) / intensities.size();
        }

        class MSstatsLine
        {
        public :
            MSstatsLine(
                    bool _has_fraction,
                    const String& _accession,
                    const String& _sequence,
                    const String& _precursor_charge,
                    const String& _fragment_ion,
                    const String& _frag_charge,
                    const String& _isotope_label_type,
                    const String& _condition,
                    const String& _bioreplicate,
                    const String& _run,
                    const String& _fraction
            ): has_fraction_(_has_fraction),
               accession_(_accession),
               sequence_(_sequence),
               precursor_charge_(_precursor_charge),
               fragment_ion_(_fragment_ion),
               frag_charge_(_frag_charge),
               isotope_label_type_(_isotope_label_type),
               condition_(_condition),
               bioreplicate_(_bioreplicate),
               run_(_run),
               fraction_(_fraction) {}
            
            const String& accession() const {return this->accession_;}
            const String& sequence() const {return this->sequence_;}
            const String& precursor_charge() const {return this->precursor_charge_;}
            const String& run() const {return this->run_;}

            String toString() const
            {
              const String delim(",");
              return  accession_
                      + delim + sequence_
                      + delim + precursor_charge_
                      + delim + fragment_ion_
                      + delim + frag_charge_
                      + delim + isotope_label_type_
                      + delim + condition_
                      + delim + bioreplicate_
                      + delim + run_
                      + (this->has_fraction_ ? delim + String(fraction_) : "");
            }

            friend bool operator<(const MSstatsLine &l,
                                  const MSstatsLine &r) {

              return std::tie(l.accession_, l.run_, l.condition_, l.bioreplicate_, l.precursor_charge_, l.sequence_) <
                     std::tie(r.accession_, r.run_, r.condition_, r.bioreplicate_, r.precursor_charge_, r.sequence_);
            }


        private:
            bool has_fraction_;
            String accession_;
            String sequence_;
            String precursor_charge_;
            String fragment_ion_;
            String frag_charge_;
            String isotope_label_type_;
            String condition_;
            String bioreplicate_;
            String run_;
            String fraction_;
        };
 
        class MSstatsTMTLine
        {
        public :
            MSstatsTMTLine(
                    bool _has_fraction,
                    const String& _accession,
                    const String& _sequence,
                    const String& _precursor_charge,
                    const String& _channel,
                    const String& _condition,
                    const String& _bioreplicate,
                    const String& _run,
                    const String& _mixture,
                    const String& _techmixture,
                    const String& _fraction
            ): has_fraction_(_has_fraction),
               accession_(_accession),
               sequence_(_sequence),
               precursor_charge_(_precursor_charge),
               channel_(_channel),
               condition_(_condition),
               bioreplicate_(_bioreplicate),
               run_(_run),
               mixture_(_mixture),
               techmixture_(_techmixture),
               fraction_(_fraction) {}

            const String& accession() const {return this->accession_;}
            const String& sequence() const {return this->sequence_;}
            const String& precursor_charge() const {return this->precursor_charge_;}
            const String& run() const {return this->run_;}

            String toString() const
            {
              const String delim(",");
              return  accession_
                      + delim + sequence_
                      + delim + precursor_charge_
                      + delim + channel_
                      + delim + condition_
                      + delim + bioreplicate_
                      + delim + run_
                      + delim + mixture_
                      + delim + techmixture_
                      + (this->has_fraction_ ? delim + String(fraction_) : "");
            }

            friend bool operator<(const MSstatsTMTLine &l,
                                  const MSstatsTMTLine &r) {

              return std::tie(l.accession_, l.run_, l.condition_, l.bioreplicate_, l.mixture_, l.precursor_charge_, l.sequence_) <
                     std::tie(r.accession_, r.run_, r.condition_, r.bioreplicate_, r.mixture_, r.precursor_charge_, r.sequence_);
            }


        private:
            bool has_fraction_;
            String accession_;
            String sequence_;
            String precursor_charge_;
            String channel_;
            String condition_;
            String bioreplicate_;
            String run_;
            String mixture_;
            String techmixture_;
            String fraction_;
        };
        
     }; 
} // namespace OpenMS
