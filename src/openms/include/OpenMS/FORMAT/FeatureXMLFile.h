// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_FEATUREXMLFILE_H
#define OPENMS_FORMAT_FEATUREXMLFILE_H

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <iosfwd>

namespace OpenMS
{
  /**
  @brief This class provides Input/Output functionality for feature maps

      A documented schema for this format can be found at http://open-ms.sourceforge.net/schemas/.

  @todo Take care that unique ids are assigned properly by TOPP tools before calling FeatureXMLFile::store().  There will be a message on LOG_INFO but we will make no attempt to fix the problem in this class.  (all developers)

  @note This format will eventually be replaced by the HUPO-PSI AnalysisXML (mzIdentML and mzQuantML) formats!

  @ingroup FileIO
*/
  class OPENMS_DLLAPI FeatureXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile,
    public ProgressLogger
  {

public:

    /** @name Constructors and Destructor */
    //@{
    ///Default constructor
    FeatureXMLFile();
    ///Destructor
    ~FeatureXMLFile();
    //@}

    /**
        @brief loads the file with name @p filename into @p map and calls updateRanges().

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String & filename, FeatureMap & feature_map);

    Size loadSize(const String & filename);

    /**
        @brief stores the map @p feature_map in file with name @p filename.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String & filename, const FeatureMap & feature_map);

    /// Mutable access to the options for loading/storing
    FeatureFileOptions & getOptions();

    /// Non-mutable access to the options for loading/storing
    const FeatureFileOptions & getOptions() const;

    /// setter for options for loading/storing
    void setOptions(const FeatureFileOptions &);

protected:

    // restore default state for next load/store operation
    void resetMembers_();

    // Docu in base class
    virtual void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname);

    // Docu in base class
    virtual void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes);

    // Docu in base class
    virtual void characters(const XMLCh * const chars, const XMLSize_t length);

    /// Writes a feature to a stream
    void writeFeature_(const String & filename, std::ostream & os, const Feature & feat, const String & identifier_prefix, UInt64 identifier, UInt indentation_level);

    /// Writes a peptide identification to a stream (for assigned/unassigned peptide identifications)
    void writePeptideIdentification_(const String & filename, std::ostream & os, const PeptideIdentification & id, const String & tag_name, UInt indentation_level);


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
    Feature * current_feature_;
    /// Feature map pointer for reading
    FeatureMap* map_;
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
    MetaInfoInterface * last_meta_;

    /// Temporary protein ProteinIdentification
    ProteinIdentification prot_id_;
    /// Temporary peptide ProteinIdentification
    PeptideIdentification pep_id_;
    /// Temporary protein hit
    ProteinHit prot_hit_;
    /// Temporary peptide hit
    PeptideHit pep_hit_;
    /// Map from protein id to accession
    Map<String, String> proteinid_to_accession_;
    /// Map from search identifier concatenated with protein accession to id
    Map<String, Size> accession_to_id_;
    /// Map from identification run identifier to file xs:id (for linking peptide identifications to the corresponding run)
    Map<String, String> identifier_id_;
    /// Map from file xs:id to identification run identifier (for linking peptide identifications to the corresponding run)
    Map<String, String> id_identifier_;
    /// Temporary search parameters file
    ProteinIdentification::SearchParameters search_param_;

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_FeatureXMLFile_H
