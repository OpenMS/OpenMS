// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer, Axel Walter $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <vector>

namespace OpenMS
{
  class ConsensusMap;
  class FeatureMap;

  /**
      @brief File adapter for QcML files used to load and store QcML files

      This Class is supposed to internally collect the data for the qcML File

      A documented schema for this format can be found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

      @ingroup FileIO
  */
  class OPENMS_DLLAPI QcMLFile :
    public Internal::XMLHandler,
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    /// Representation of a quality parameter
    class OPENMS_DLLAPI QualityParameter
    {
    public:
      std::string name; ///< Name
      std::string id; ///< Identifier
      std::string value; ///< Value
      std::string cvRef; ///< cv reference
      std::string cvAcc; ///< cv accession
      std::string unitRef; ///< cv reference of the unit
      std::string unitAcc; ///< cv accession of the unit
      std::string flag; ///< cv accession of the unit

      ///Default constructor
      QualityParameter();

      QualityParameter(const QualityParameter& rhs);

      QualityParameter& operator=(const QualityParameter& rhs);
      bool operator==(const QualityParameter& rhs) const;
      bool operator<(const QualityParameter& rhs) const;
      bool operator>(const QualityParameter& rhs) const;

      std::string toXMLString(UInt indentation_level) const;
    };

    /// Representation of an attachment
    class OPENMS_DLLAPI Attachment
    {
    public:
      std::string name; ///< Name
      std::string id; ///< Name
      std::string value; ///< Value
      std::string cvRef; ///< cv reference
      std::string cvAcc; ///< cv accession
      std::string unitRef; ///< cv reference of the unit
      std::string unitAcc; ///< cv accession of the unit
      std::string binary; ///< binary content of the attachment
      std::string qualityRef; ///< reference to qp to which attachment, if empty attached to run/set
      std::vector<std::string> colTypes; ///< type of the cols if QP has a table of values
      std::vector< std::vector<std::string> > tableRows; ///< cell values if QP has a table, type see colType
      //~ TODO -schema- coltypes with full definition (uintRef, unitAcc)

      ///Default constructor
      Attachment();

      Attachment(const Attachment& rhs);

      Attachment& operator=(const Attachment& rhs);
      bool operator==(const Attachment& rhs) const;
      bool operator<(const Attachment& rhs) const;
      bool operator>(const Attachment& rhs) const;

      std::string toXMLString(UInt indentation_level) const;
      std::string toCSVString(const std::string& separator) const;
    };

    ///Default constructor
    QcMLFile();
    ///Destructor
    ~QcMLFile() override;

    std::string map2csv(const std::map< std::string, std::map<std::string, std::string> >& cvs_table, const std::string& separator) const;
    std::string exportIDstats(const std::string& filename) const;
    
    /// Registers a run in the qcml file with the respective mappings
    void registerRun(const std::string& id, const std::string& name);
    /// Registers a set in the qcml file with the respective mappings
    void registerSet(const std::string& id, const std::string& name, const std::set<std::string>& names);
    /// Just adds a qualityparameter to run by the name r
    void addRunQualityParameter(const std::string& r, const QualityParameter& qp);
    /// Just adds a attachment to run by the name r
    void addRunAttachment(const std::string& r, const Attachment& at);
    /// Just adds a qualityparameter to set by the name r
    void addSetQualityParameter(const std::string& r, const QualityParameter& qp);
    /// Just adds a attachment to set by the name r
    void addSetAttachment(const std::string& r, const Attachment& at);
    /// Removes attachments referencing a id given in ids, from run/set r. All attachments if no attachment name is given with at.
    void removeAttachment(const std::string& r, std::vector<std::string>& ids, const std::string& at = "");
    /// Removes attachment with cv accession at from run/set r.
    void removeAttachment(const std::string& r, const std::string& at);
    /// Removes attachment with cv accession at from  all runs/sets.
    void removeAllAttachments(const std::string& at);
    /// Just removes qualityparameter going by one of the ID attributes given in ids.
    void removeQualityParameter(const std::string& r, std::vector<std::string>& ids);
    /// merges the given QCFile into this one
    void merge(const QcMLFile & addendum, const std::string& setname = "");
    /// collects the values of given QPs (as CVid) of the given set
    void/* std::vector<std::string>& */ collectSetParameter(const std::string& setname, const std::string& qp, std::vector<std::string>& ret);
    /// Returns a std::string of a tab separated rows if found empty string else from run/set by the name filename of the qualityparameter by the name qpname
    std::string exportAttachment(const std::string& filename, const std::string& qpname) const; 
    /// Returns a std::string value in quotation of a qualityparameter by the name qpname in run/set by the name filename
    std::string exportQP(const std::string& filename, const std::string& qpname) const;
    /// Returns a std::string of a tab separated qualityparameter by the name qpname in run/set by the name filename
    std::string exportQPs(const std::string& filename, const StringList& qpnames) const;
    /// Gives the ids of the registered runs in the vector ids.
    void getRunIDs (std::vector<std::string>& ids) const;
    /// Gives the names of the registered runs in the vector ids.
    void getRunNames (std::vector<std::string>& ids) const;
    /// Returns true if the given run id is present in this file, if checkname is true it also checks the names
    bool existsRun(const std::string& filename, bool checkname = false) const;
    /// Returns true if the given set id is present in this file, if checkname is true it also checks the names
    bool existsSet(const std::string& filename, bool checkname = false) const;
    /// Returns the ids of the parameter name given if found in given run empty else
    void existsRunQualityParameter(const std::string& filename, const std::string& qpname, std::vector<std::string>& ids) const;
    /// Returns the ids of the parameter name given if found in given set, empty else
    void existsSetQualityParameter(const std::string& filename, const std::string& qpname, std::vector<std::string>& ids) const;
    /// Calculation and collection of QC data
    /**
      @brief Collects QC data in qualityParameters and qualityAttachments
      @param[in] prot_ids protein identifications from ID file
      @param[in] pep_ids peptide identifications
      @param[in] feature_map FeatureMap from feature file (featureXML)
      @param[in] consensus_map ConsensusMap from consensus file (consensusXML)
      @param[in] inputfile_raw mzML input file name
      @param[in] remove_duplicate_features removes duplicates in a set of merged features
      @param[in] exp MSExperiment to extract QC data from, prior sortSpectra() and updateRanges() required
    */
    void collectQCData(std::vector<ProteinIdentification>& prot_ids,
                       PeptideIdentificationList& pep_ids,
                       const FeatureMap& feature_map,
                       const ConsensusMap& consensus_map,
                       const std::string& inputfile_raw,
                       const bool remove_duplicate_features,
                       const MSExperiment& exp);
    ///Store the QCFile
    /**
      @brief Store the qcML file
      @param[out] filename qcML output file name
    */
    void store(const std::string& filename) const;

    ///Load a QCFile
    void load(const std::string & filename);

    //~ int siz; //debug

protected:
    // Docu in base class
    void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname) override;

    // Docu in base class
    void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

    // Docu in base class
    void characters(const XMLCh * const chars, const XMLSize_t length) override;

    std::map<std::string, std::vector< QualityParameter > > runQualityQPs_; //TODO run name attribute to schema of RunQuality
    std::map<std::string, std::vector< Attachment > > runQualityAts_;
    std::map<std::string, std::vector< QualityParameter > > setQualityQPs_;
    std::map<std::string, std::vector< Attachment > > setQualityAts_;
    std::map<std::string, std::set< std::string > > setQualityQPs_members_;
    std::map<std::string, std::string > run_Name_ID_map_;
    std::map<std::string, std::string > set_Name_ID_map_;

    std::string tag_;
    UInt progress_;
    QualityParameter qp_;
    Attachment at_;
    std::vector<std::string> row_;
    std::vector<std::string> header_;
    std::string name_;
    std::string run_id_;
    std::set<std::string> names_;
    std::vector<QualityParameter> qps_;
    std::vector<Attachment> ats_;
  };

} // namespace OpenMS
