// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <iosfwd>
#include <map>

namespace OpenMS
{
  class Feature;
  class FeatureMap;

  namespace Internal
  {

  /**
    @brief This class provides Input/Output functionality for feature maps

    A documented schema for this format can be found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

    @todo Take care that unique ids are assigned properly by TOPP tools before
    calling FeatureXMLFile::store().  There will be a message on OPENMS_LOG_INFO but
    we will make no attempt to fix the problem in this class.  (all developers)

    @ingroup FileIO
  */
  class OPENMS_DLLAPI FeatureXMLHandler :
    public Internal::XMLHandler,
    public ProgressLogger
  {

public:

    /** @name Constructors and Destructor */
    //@{
    ///Default constructor
    FeatureXMLHandler(FeatureMap& map, const String& filename);
    FeatureXMLHandler(const FeatureMap& map, const String& filename);
    ///Destructor
    ~FeatureXMLHandler() override;
    //@}


    /// Docu in base class XMLHandler::writeTo
    void writeTo(std::ostream& os) override;

    /// Mutable access to the options for loading/storing
    FeatureFileOptions& getOptions();

    /// Non-mutable access to the options for loading/storing
    const FeatureFileOptions& getOptions() const;

    /// setter for options for loading/storing
    void setOptions(const FeatureFileOptions&);

    void setSizeOnly(const bool size_only)
    {
      size_only_ = size_only;
    }

    Size getSize() const
    {
      return expected_size_;
    }

protected:

    // restore default state for next load/store operation
    void resetMembers_();

    // Docu in base class
    void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

    // Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

    // Docu in base class
    void characters(const XMLCh* const chars, const XMLSize_t length) override;

    /// Writes a feature to a stream
    void writeFeature_(const String& filename, std::ostream& os, const Feature& feat, const String& identifier_prefix, UInt64 identifier, UInt indentation_level);

    /// Writes a peptide identification to a stream (for assigned/unassigned peptide identifications)
    void writePeptideIdentification_(const String& filename, std::ostream& os, const PeptideIdentification& id, const String& tag_name, UInt indentation_level);

    /**
        @brief update the pointer to the current feature

        @param create If true, a new (empty) Feature is added at the appropriate subordinate_feature_level_
    */
    void updateCurrentFeature_(bool create);

    /// allows for early return in parsing functions when certain sections should be ignored
    /// <=0 - parsing ON
    ///  >0 - this number of tags have been entered that forbid parsing and need to be exited before parsing continues
    Int disable_parsing_;

    /// points to the last open &lt;feature&gt; tag (possibly a subordinate feature)
    Feature* current_feature_;
    /// Feature map pointer for writing
    FeatureMap* map_;
    /// Feature map pointer for reading
    const FeatureMap* cmap_;
    /// Options that can be set
    FeatureFileOptions options_;
    /// only parse until "count" tag is reached (used in loadSize())
    bool size_only_;
    /// holds the putative size given in count
    Size expected_size_;

    /**@name temporary data structures to hold parsed data */
    //@{
    Param param_;
    ConvexHull2D::PointArrayType current_chull_;
    DPosition<2> hull_position_;
    //@}

    /// current dimension of the feature position, quality, or convex hull point
    UInt dim_;

    /// for downward compatibility, all tags in the old description must be ignored
    bool in_description_;

    /// level in Feature stack during parsing
    Int subordinate_feature_level_;

    /// Pointer to last read object as a MetaInfoInterface, or null.
    MetaInfoInterface* last_meta_;

    /// Temporary protein ProteinIdentification
    ProteinIdentification prot_id_;
    /// Temporary peptide ProteinIdentification
    PeptideIdentification pep_id_;
    /// Temporary protein hit
    ProteinHit prot_hit_;
    /// Temporary peptide hit
    PeptideHit pep_hit_;
    /// Map from protein id to accession
    std::map<String, String> proteinid_to_accession_;
    /// Map from search identifier concatenated with protein accession to id
    std::map<String, Size> accession_to_id_;
    /// Map from identification run identifier to file xs:id (for linking peptide identifications to the corresponding run)
    std::map<String, String> identifier_id_;
    /// Map from file xs:id to identification run identifier (for linking peptide identifications to the corresponding run)
    std::map<String, String> id_identifier_;
    /// Temporary search parameters file
    ProteinIdentification::SearchParameters search_param_;

  };

} // namespace Internal
} // namespace OpenMS

