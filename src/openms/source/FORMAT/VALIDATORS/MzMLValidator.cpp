// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>

using namespace xercesc;
using namespace std;

namespace OpenMS::Internal
{

    MzMLValidator::MzMLValidator(const CVMappings & mapping, const ControlledVocabulary & cv) :
      SemanticValidator(mapping, cv),
      binary_data_array_(),
      binary_data_type_()
    {
      setCheckUnits(true);
    }

    MzMLValidator::~MzMLValidator()
    = default;

    //This method needed to be reimplemented to
    // - check CV term values
    // - handle referenceableParamGroups
    // - check if binaryDataArray name and type match
    void MzMLValidator::startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const Attributes & attributes)
    {
      String tag = sm_.convert(qname);
      String parent_tag;
      if (!open_tags_.empty())
      {
        parent_tag = open_tags_.back();
      }
      String path = getPath_() + "/" + cv_tag_ + "/@" + accession_att_;
      open_tags_.push_back(tag);

      if (tag == "referenceableParamGroup")
      {
        current_id_ = attributeAsString_(attributes, "id");
      }
      else if (tag == "referenceableParamGroupRef")
      {
        const std::vector<CVTerm> & terms = param_groups_[attributeAsString_(attributes, "ref")];
        for (Size i = 0; i < terms.size(); ++i)
        {
          handleTerm_(path, terms[i]);
        }
      }
      else if (tag == "binaryDataArray")
      {
        binary_data_array_ = "";
        binary_data_type_ = "";
      }
      else if (tag == cv_tag_)
      {
        //extract accession, name and value
        CVTerm parsed_term;
        getCVTerm_(attributes, parsed_term);

        //check if the term is unknown
        if (!cv_.exists(parsed_term.accession))
        {
          warnings_.push_back(String("Unknown CV term: '") + parsed_term.accession + " - " + parsed_term.name + "' at element '" + getPath_(1) + "'");
          return;
        }

        //check if the term is obsolete
        if (cv_.getTerm(parsed_term.accession).obsolete)
        {
          warnings_.push_back(String("Obsolete CV term: '") + parsed_term.accession + " - " + parsed_term.name + "' at element '" + getPath_(1) + "'");
        }

        //actual handling of the term
        if (parent_tag == "referenceableParamGroup")
        {
          param_groups_[current_id_].push_back(parsed_term);
        }
        else
        {
          handleTerm_(path, parsed_term);
        }
      }
    }

    //reimplemented in order to remove the "indexedmzML" tag from the front (if present)
    String MzMLValidator::getPath_(UInt remove_from_end) const
    {
      String path;
      if (!open_tags_.empty() && open_tags_.front() == "indexedmzML")
      {
        path.concatenate(open_tags_.begin() + 1, open_tags_.end() - remove_from_end, "/");
      }
      else
      {
        path.concatenate(open_tags_.begin(), open_tags_.end() - remove_from_end, "/");
      }
      path = String("/") + path;
      return path;
    }

    //reimplemented to
    // - catch non-PSI CVs
    // - check if binaryDataArray name and type match
    void MzMLValidator::handleTerm_(const String & path, const CVTerm & parsed_term)
    {
      //some CVs cannot be validated because they use 'part_of' which spoils the inheritance
      if (parsed_term.accession.hasPrefix("GO:"))
      {
        return;
      }
      if (parsed_term.accession.hasPrefix("BTO:"))
      {
        return;
      }
      //check binary data array terms
      if (path.hasSuffix("/binaryDataArray/cvParam/@accession"))
      {
        //binary data array
        if (cv_.isChildOf(parsed_term.accession, "MS:1000513"))
        {
          binary_data_array_ = parsed_term.accession;
        }
        //binary data type
        if (cv_.isChildOf(parsed_term.accession, "MS:1000518"))
        {
          binary_data_type_ = parsed_term.accession;
        }
        //if both are parsed, check if they match
        if (!binary_data_type_.empty() && !binary_data_array_.empty())
        {
          if (!ListUtils::contains(cv_.getTerm(binary_data_array_).xref_binary, binary_data_type_))
          {
            errors_.push_back(String("Binary data array of type '") + binary_data_array_ + " ! " + cv_.getTerm(binary_data_array_).name + "' cannot have the value type '" + binary_data_type_ + " ! " + cv_.getTerm(binary_data_type_).name + "'.");
          }
        }
      }

      SemanticValidator::handleTerm_(path, parsed_term);
    }

} // namespace OpenMS   // namespace Internal
