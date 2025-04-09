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
      String name; ///< Name
      String id; ///< Identifier
      String value; ///< Value
      String cvRef; ///< cv reference
      String cvAcc; ///< cv accession
      String unitRef; ///< cv reference of the unit
      String unitAcc; ///< cv accession of the unit
      String flag; ///< cv accession of the unit

      ///Default constructor
      QualityParameter();

      QualityParameter(const QualityParameter& rhs);

      QualityParameter& operator=(const QualityParameter& rhs);
      bool operator==(const QualityParameter& rhs) const;
      bool operator<(const QualityParameter& rhs) const;
      bool operator>(const QualityParameter& rhs) const;

      String toXMLString(UInt indentation_level) const;
    };

    /// Representation of an attachment
    class OPENMS_DLLAPI Attachment
    {
    public:
      String name; ///< Name
      String id; ///< Name
      String value; ///< Value
      String cvRef; ///< cv reference
      String cvAcc; ///< cv accession
      String unitRef; ///< cv reference of the unit
      String unitAcc; ///< cv accession of the unit
      String binary; ///< binary content of the attachment
      String qualityRef; ///< reference to qp to which attachment, if empty attached to run/set
      std::vector<String> colTypes; ///< type of the cols if QP has a table of values
      std::vector< std::vector<String> > tableRows; ///< cell values if QP has a table, type see colType
      //~ TODO -schema- coltypes with full definition (uintRef, unitAcc)

      ///Default constructor
      Attachment();

      Attachment(const Attachment& rhs);

      Attachment& operator=(const Attachment& rhs);
      bool operator==(const Attachment& rhs) const;
      bool operator<(const Attachment& rhs) const;
      bool operator>(const Attachment& rhs) const;

      String toXMLString(UInt indentation_level) const;
      String toCSVString(const String& separator) const;
    };

    ///Default constructor
    QcMLFile();
    ///Destructor
    ~QcMLFile() override;

    String map2csv(const std::map< String, std::map<String, String> >& cvs_table, const String& separator) const;
    String exportIDstats(const String& filename) const;
    
    /// Registers a run in the qcml file with the respective mappings
    void registerRun(const String& id, const String& name);
    /// Registers a set in the qcml file with the respective mappings
    void registerSet(const String& id, const String& name, const std::set<String>& names);
    /// Just adds a qualityparameter to run by the name r
    void addRunQualityParameter(const String& r, const QualityParameter& qp);
    /// Just adds a attachment to run by the name r
    void addRunAttachment(const String& r, const Attachment& at);
    /// Just adds a qualityparameter to set by the name r
    void addSetQualityParameter(const String& r, const QualityParameter& qp);
    /// Just adds a attachment to set by the name r
    void addSetAttachment(const String& r, const Attachment& at);
    /// Removes attachments referencing a id given in ids, from run/set r. All attachments if no attachment name is given with at.
    void removeAttachment(const String& r, std::vector<String>& ids, const String& at = "");
    /// Removes attachment with cv accession at from run/set r.
    void removeAttachment(const String& r, const String& at);
    /// Removes attachment with cv accession at from  all runs/sets.
    void removeAllAttachments(const String& at);
    /// Just removes qualityparameter going by one of the ID attributes given in ids.
    void removeQualityParameter(const String& r, std::vector<String>& ids);
    /// merges the given QCFile into this one
    void merge(const QcMLFile & addendum, const String& setname = "");
    /// collects the values of given QPs (as CVid) of the given set
    void/* std::vector<String>& */ collectSetParameter(const String& setname, const String& qp, std::vector<String>& ret);
    /// Returns a String of a tab separated rows if found empty string else from run/set by the name filename of the qualityparameter by the name qpname
    String exportAttachment(const String& filename, const String& qpname) const; 
    /// Returns a String value in quotation of a qualityparameter by the name qpname in run/set by the name filename
    String exportQP(const String& filename, const String& qpname) const;
    /// Returns a String of a tab separated qualityparameter by the name qpname in run/set by the name filename
    String exportQPs(const String& filename, const StringList& qpnames) const;
    /// Gives the ids of the registered runs in the vector ids.
    void getRunIDs (std::vector<String>& ids) const;
    /// Gives the names of the registered runs in the vector ids.
    void getRunNames (std::vector<String>& ids) const;
    /// Returns true if the given run id is present in this file, if checkname is true it also checks the names
    bool existsRun(const String& filename, bool checkname = false) const;
    /// Returns true if the given set id is present in this file, if checkname is true it also checks the names
    bool existsSet(const String& filename, bool checkname = false) const;
    /// Returns the ids of the parameter name given if found in given run empty else
    void existsRunQualityParameter(const String& filename, const String& qpname, std::vector<String>& ids) const;
    /// Returns the ids of the parameter name given if found in given set, empty else
    void existsSetQualityParameter(const String& filename, const String& qpname, std::vector<String>& ids) const;
    /// Calculation and collection of QC data
    /**
      @brief Collects QC data in qualityParameters and qualityAttachments
      @param prot_ids protein identifications from ID file
      @param pep_ids peptide identifications
      @param feature_map FeatureMap from feature file (featureXML)
      @param consensus_map ConsensusMap from consensus file (consensusXML)
      @param inputfile_raw mzML input file name
      @param remove_duplicate_features removes duplicates in a set of merged features
      @param exp MSExperiment to extract QC data from, prior sortSpectra() and updateRanges() required
    */
    void collectQCData(std::vector<ProteinIdentification>& prot_ids,
                       std::vector<PeptideIdentification>& pep_ids,
                       const FeatureMap& feature_map,
                       const ConsensusMap& consensus_map,
                       const String& inputfile_raw,
                       const bool remove_duplicate_features,
                       const MSExperiment& exp);
    ///Store the QCFile
    /**
      @brief Store the qcML file
      @param filename qcML output file name
    */
    void store(const String& filename) const;

    ///Load a QCFile
    void load(const String & filename);

    //~ int siz; //debug

protected:
    // Docu in base class
    void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname) override;

    // Docu in base class
    void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

    // Docu in base class
    void characters(const XMLCh * const chars, const XMLSize_t length) override;

    std::map<String, std::vector< QualityParameter > > runQualityQPs_; //TODO run name attribute to schema of RunQuality
    std::map<String, std::vector< Attachment > > runQualityAts_;
    std::map<String, std::vector< QualityParameter > > setQualityQPs_;
    std::map<String, std::vector< Attachment > > setQualityAts_;
    std::map<String, std::set< String > > setQualityQPs_members_;
    std::map<String, String > run_Name_ID_map_;
    std::map<String, String > set_Name_ID_map_;

    String tag_;
    UInt progress_;
    QualityParameter qp_;
    Attachment at_;
    std::vector<String> row_;
    std::vector<String> header_;
    String name_;
    String run_id_;
    std::set<String> names_;
    std::vector<QualityParameter> qps_;
    std::vector<Attachment> ats_;
  };

} // namespace OpenMS
