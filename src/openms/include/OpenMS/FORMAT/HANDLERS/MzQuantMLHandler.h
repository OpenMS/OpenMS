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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MZQUANTMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZQUANTMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/METADATA/MSQuantifications.h>

namespace OpenMS
{
  class ProgressLogger;

  namespace Internal
  {

    /**
        @brief XML handler for MzQuantMLFile

        @note Do not use this class. It is only needed in MzQuantMLFile.
    */
    class OPENMS_DLLAPI MzQuantMLHandler :
      public XMLHandler
    {
public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      MzQuantMLHandler(const MSQuantifications & msq, const String & filename, const String & version, const ProgressLogger & logger);

      /// Constructor for a read-only handler
      MzQuantMLHandler(MSQuantifications & msq, const String & filename, const String & version, const ProgressLogger & logger);

      /// Destructor
      ~MzQuantMLHandler() override;
      //@}


      // Docu in base class
      void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname) override;

      // Docu in base class
      void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

      // Docu in base class
      void characters(const XMLCh * const chars, const XMLSize_t length) override;

      //Docu in base class
      void writeTo(std::ostream & os) override;

protected:

      /// Progress logger
      const ProgressLogger & logger_;

      /// Controlled vocabulary (hopefully the psi-pi from OpenMS/share/OpenMS/CV/psi-pi.obo)
      ControlledVocabulary cv_;

      String tag_;

      MSQuantifications * msq_;

      const MSQuantifications * cmsq_;

      /// Handles CV terms
      void handleCVParam_(const String & parent_parent_tag, const String & parent_tag, const String & accession, const String & name, const String & value, const xercesc::Attributes & attributes, const String & cv_ref, const String & unit_accession = "");

      /// Handles user terms
      void handleUserParam_(const String & parent_parent_tag, const String & parent_tag, const String & name, const String & type, const String & value);

      /// Write CV term
      //~ TODO rewirte writeCVParams_ in baseclass to be more convenient
      //~ void MzQuantMLHandler::writeCVParams_(std::ostream& os, const Map< String, std::vector < CVTerm > > & , UInt indent);
      void writeCVParams_(String & s, const Map<String, std::vector<CVTerm> > &, UInt indent);

      /// Writes user terms
      void writeUserParams_(std::ostream & os, const MetaInfoInterface & meta, UInt indent);
      void writeUserParams_(String & s, const MetaInfoInterface & meta, UInt indent);

      /// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
      ControlledVocabulary::CVTerm getChildWithName_(const String & parent_accession, const String & name) const;


      /// Helper method that writes the featuremaps
      void writeFeature_(String & feature_xml, const std::vector<FeatureMap >& fm, UInt indentation_level);

      /// Helper method that writes a source file
      //void writeSourceFile_(std::ostream& os, const String& id, const SourceFile& software);


private:
      MzQuantMLHandler();
      MzQuantMLHandler(const MzQuantMLHandler & rhs);
      MzQuantMLHandler & operator=(const MzQuantMLHandler & rhs);

      std::map<String, std::vector<ExperimentalSettings> > current_files_;           // 1.rawfilesgroup_ref 2.inputfiles for each assay as ExperimentalSettings
      String current_id_;
      String current_cf_id_;
      Size current_count_;

      std::vector<MetaInfo> up_stack_;
      std::vector<CVTerm> cvp_stack_;
      MSQuantifications::Assay current_assay_;

      std::multimap<String, String> cm_cf_ids_;
      std::map<String, String> f_cf_ids_;
      std::map<String, ConsensusFeature> cf_cf_obj_;
      std::map<String, FeatureHandle> f_f_obj_;
      std::map<String, ConsensusFeature::Ratio> r_rtemp_;
      std::map<String, String> numden_r_ids_;
      std::map<String, ConsensusFeature::Ratio> r_r_obj_;

      //~ Software current_sw_;
      std::map<String, Software> current_sws_;
      std::map<int, DataProcessing> current_orderedps_;
      std::pair<int, DataProcessing> current_dp_;
      std::set<DataProcessing::ProcessingAction> current_pas_;

      std::vector<String> current_col_types_;
      std::vector<double> current_dm_values_;
      std::vector<double> current_row_;

    };
  }   // namespace Internal
} // namespace OpenMS

#endif
