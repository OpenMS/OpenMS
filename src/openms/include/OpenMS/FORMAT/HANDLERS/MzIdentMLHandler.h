// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <vector>
#include <map>
#include <set>

namespace OpenMS
{
  class ProgressLogger;

  namespace Internal
  {

  /**
    @brief Represents a object which can store the information of an analysisXML instance

    @ingroup Metadata
  */
  class OPENMS_DLLAPI IdentificationHit :
    public MetaInfoInterface
  {
  public:
    /// @name Constructors, Destructors, Assignment Operators
    //@{
    /// Default constructor
    IdentificationHit() = default;
    
    /// Copy constructor
    IdentificationHit(const IdentificationHit&) = default;
    
    /// Virtual destructor
    virtual ~IdentificationHit() = default;
    
    /// Move constructor
    IdentificationHit(IdentificationHit&&) noexcept = default;
    
    /// Copy assignment operator
    IdentificationHit& operator=(const IdentificationHit&) = default;
    
    /// Move assignment operator
    IdentificationHit& operator=(IdentificationHit&&) noexcept = default;
    //@}

    /// @name Equality and Inequality Operators
    //@{
    /// Checks for equality with another IdentificationHit object
    bool operator==(const IdentificationHit& rhs) const noexcept;
    
    /// Checks for inequality with another IdentificationHit object
    bool operator!=(const IdentificationHit& rhs) const noexcept;
    //@}

    /// @name Accessors
    //@{
    /// Sets the identifier
    void setId(const std::string& id) noexcept;
    
    /// Returns the identifier
    const std::string& getId() const noexcept;
    
    /// Sets the charge state of the peptide
    void setCharge(int charge) noexcept;
    
    /// Returns the charge state of the peptide
    int getCharge() const noexcept;
    
    /// Sets the calculated mass to charge ratio
    void setCalculatedMassToCharge(double mz) noexcept;
    
    /// Returns the calculated mass to charge ratio
    double getCalculatedMassToCharge() const noexcept;
    
    /// Sets the experimental mass to charge ratio
    void setExperimentalMassToCharge(double mz) noexcept;
    
    /// Returns the experimental mass to charge ratio
    double getExperimentalMassToCharge() const noexcept;
    
    /// Sets the name
    void setName(const std::string& name) noexcept;
    
    /// Returns the name
    const std::string& getName() const noexcept;
    
    /// Sets whether the peptide passed the threshold
    void setPassThreshold(bool pass) noexcept;
    
    /// Returns whether the peptide passed the threshold
    bool getPassThreshold() const noexcept;
    
    /// Sets the rank of the peptide
    void setRank(int rank) noexcept;
    
    /// Returns the rank of the peptide
    int getRank() const noexcept;
    //@}

  private:
    std::string id_;                              ///< Identifier
    int charge_ = 0;                             ///< Peptide charge
    double calculated_mass_to_charge_ = 0.0;     ///< Calculated mass to charge ratio
    double experimental_mass_to_charge_ = 0.0;   ///< Experimental mass to charge ratio
    std::string name_;                           ///< Name
    bool pass_threshold_ = true;                 ///< Pass threshold
    int rank_ = 0;                               ///< Rank of the peptide
  };

  /**
    @brief Represents a object which can store the information of an analysisXML instance

        //@todo docu (Andreas)

        @ingroup Metadata
  */
  class OPENMS_DLLAPI SpectrumIdentification :
    public MetaInfoInterface
  {
public:

    /// @name constructors,destructors,assignment operator
    //@{
    /// Default constructor
    SpectrumIdentification() = default;
    /// Destructor
    virtual ~SpectrumIdentification();
    /// Copy constructor
    SpectrumIdentification(const SpectrumIdentification &) = default;
    /// Move constructor
    SpectrumIdentification(SpectrumIdentification&&) = default;
    /// Assignment operator
    SpectrumIdentification & operator=(const SpectrumIdentification &) = default;
    /// Move assignment operator
    SpectrumIdentification& operator=(SpectrumIdentification&&) & = default;
    /// Equality operator
    bool operator==(const SpectrumIdentification & rhs) const;
    /// Inequality operator
    bool operator!=(const SpectrumIdentification & rhs) const;
    //@}

    // @name Accessors
    //@{
    /// sets the identification hits of this spectrum identification (corresponds to single peptide hit in the list)
    void setHits(const std::vector<IdentificationHit> & hits);

    /// adds a single identification hit to the hits
    void addHit(const IdentificationHit & hit);

    /// returns the identification hits of this spectrum identification
    const std::vector<IdentificationHit> & getHits() const;
    //@}

protected:

    std::string id_; ///< Identifier
    std::vector<IdentificationHit> hits_; ///< Single peptide hits
  };

    /**
      @brief Represents a object which can store the information of an analysisXML instance

          //@todo docu (Andreas)

          @ingroup Metadata
    */
    class OPENMS_DLLAPI Identification :
      public MetaInfoInterface
    {
  public:

      /// @name constructors,destructors,assignment operator
      //@{

      /// Default constructor
      Identification() = default;
      /// Copy constructor
      Identification(const Identification & source) = default;
      /// Move constructor
      Identification(Identification&&) = default;
      /// Destructor
      virtual ~Identification();

      /// Assignment operator
      Identification & operator=(const Identification & source) = default;
      /// Move assignment operator
      Identification& operator=(Identification&&) & = default;

      /// Equality operator
      bool operator==(const Identification & rhs) const;
      /// Inequality operator
      bool operator!=(const Identification & rhs) const;
      //@}

      /// @name Accessors
      //@{
      /// sets the date and time the file was written
      void setCreationDate(const DateTime & date);

      /// returns the date and time the file was created
      const DateTime & getCreationDate() const;

      /// sets the spectrum identifications
      void setSpectrumIdentifications(const std::vector<SpectrumIdentification> & ids);

      /// adds a spectrum identification
      void addSpectrumIdentification(const SpectrumIdentification & id);

      /// returns the spectrum identifications stored
      const std::vector<SpectrumIdentification> & getSpectrumIdentifications() const;
      //@}
  protected:
      std::string id_; ///< Identifier
      DateTime creation_date_; ///< Date and time the search was performed
      std::vector<SpectrumIdentification> spectrum_identifications_;
    };

    /**
        @brief XML STREAM handler for MzIdentMLFile

        In read-mode, this class will parse an MzIdentML XML file and append the input
        identifications to the provided PeptideIdentifications and ProteinIdentifications.

        @note Do not use this class. It is only needed in MzIdentMLFile.
        @note DOM and STREAM handler for MzIdentML have the same interface for legacy id structures.
    */
    class OPENMS_DLLAPI MzIdentMLHandler :
      public XMLHandler
    {
public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler for internal identification structures
      MzIdentMLHandler(const std::vector<ProteinIdentification>& pro_id, const PeptideIdentificationList& pep_id, const std::string& filename, const std::string& version, const ProgressLogger& logger);

      /// Constructor for a read-only handler for internal identification structures
      MzIdentMLHandler(std::vector<ProteinIdentification>& pro_id, PeptideIdentificationList& pep_id, const std::string& filename, const std::string& version, const ProgressLogger& logger);

      /// Destructor
      ~MzIdentMLHandler() override;
      //@}


      // Docu in base class
      void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

      // Docu in base class
      void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

      // Docu in base class
      void characters(const XMLCh* const chars, const XMLSize_t length) override;

      //Docu in base class
      void writeTo(std::ostream& os) override;

protected:
      /// Progress logger
      const ProgressLogger& logger_;

      ///Controlled vocabulary (psi-ms from OpenMS/share/OpenMS/CV/psi-ms.obo)
      ControlledVocabulary cv_;
      ///Controlled vocabulary for modifications (unimod from OpenMS/share/OpenMS/CV/unimod.obo)
      ControlledVocabulary unimod_;

      //~ PeakMap* ms_exp_;

      ///XML tag parse element
      std::string tag_;

      ///Identification Item
      Identification* id_;
      ///internal Identification Item for proteins
      std::vector<ProteinIdentification>* pro_id_;
      ///Identification Item for peptides
      PeptideIdentificationList* pep_id_;

      const Identification* cid_;
      const std::vector<ProteinIdentification>* cpro_id_;
      const PeptideIdentificationList* cpep_id_;

      ///SpectrumIdentification Item
      SpectrumIdentification current_spectrum_id_;

      ///IdentificationHit Item
      IdentificationHit current_id_hit_;

      /// Handles CV terms
      void handleCVParam_(const std::string& parent_parent_tag, const std::string& parent_tag, const std::string& accession, /* const std::string& name, */ /* const std::string& value, */ const xercesc::Attributes& attributes, const std::string& cv_ref /* ,  const std::string& unit_accession="" */);

      /// Handles user terms
      void handleUserParam_(const std::string& parent_parent_tag, const std::string& parent_tag, const std::string& name, const std::string& type, const std::string& value);

      /// Writes user terms
      void writeMetaInfos_(std::string& s, const MetaInfoInterface& meta, UInt indent) const;

      /// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
      ControlledVocabulary::CVTerm getChildWithName_(const std::string& parent_accession, const std::string& name) const;

      /// Helper method that writes a source file
      //void writeSourceFile_(std::ostream& os, const std::string& id, const SourceFile& software);

      /// Helper method that writes the Enzymes
      void writeEnzyme_(std::string& s, const DigestionEnzymeProtein& enzy, UInt miss, UInt indent) const;

      /// Helper method that writes the modification search params (fixed or variable)
      void writeModParam_(std::string& s, const std::vector<std::string>& mod_names, bool fixed, UInt indent) const;

      /// Helper method that writes the FragmentAnnotations section of a spectrum identification
      void writeFragmentAnnotations_(std::string& s, const std::vector<PeptideHit::PeakAnnotation>& annotations, UInt indent, bool is_ppxl) const;

      /// Convenience method to remove the [] from OpenMS internal file uri representation
      std::string trimOpenMSfileURI(const std::string& file) const;

      /// Abstraction of PeptideHit loop for most PeptideHits
      void writePeptideHit(const PeptideHit& hit,
                                PeptideIdentificationList::const_iterator& it,
                                std::map<std::string, std::string>& pep_ids,
                                const std::string& cv_ns, std::set<std::string>& sen_set,
                                std::map<std::string, std::string>& sen_ids,
                                std::map<std::string, std::vector<std::string> >& pep_evis,
                                std::map<std::string, double>& pp_identifier_2_thresh,
                                std::string& sidres);

      /// Abstraction of PeptideHit loop for XL-MS data from OpenPepXL
      void writeXLMSPeptideHit(const PeptideHit& hit,
                                PeptideIdentificationList::const_iterator& it,
                                const std::string& ppxl_linkid, std::map<std::string, std::string>& pep_ids,
                                const std::string& cv_ns, std::set<std::string>& sen_set,
                                std::map<std::string, std::string>& sen_ids,
                                std::map<std::string, std::vector<std::string> >& pep_evis,
                                std::map<std::string, double>& pp_identifier_2_thresh,
                                double ppxl_crosslink_mass,
                                std::map<std::string, std::string>& ppxl_specref_2_element,
                                std::string& sid, bool alpha_peptide);

private:
      MzIdentMLHandler();

      /// Load CVs and precompute the cached CV child-term set; shared by both constructors.
      void initCvCaches_();
      MzIdentMLHandler(const MzIdentMLHandler& rhs);
      MzIdentMLHandler& operator=(const MzIdentMLHandler& rhs);
      std::map<std::string, AASequence> pep_sequences_;
      std::map<std::string, std::string> pp_identifier_2_sil_; ///< mapping peptide/proteinidentification identifier_ to spectrumidentificationlist
      std::map<std::string, std::string> sil_2_sdb_; ///< mapping spectrumidentificationlist to the search data bases
      std::map<std::string, std::string> sil_2_sdat_; ///< mapping spectrumidentificationlist to the search input
      std::map<std::string, std::string> ph_2_sdat_; ///< mapping identification runs (mapping PeptideIdentifications and ProteinIdentifications via .getIdentifier()) to spectra data
      std::map<std::string, std::string> sil_2_sip_; ///< mapping spectrumidentificationlist to the search protocol (where the params are at)
      AASequence actual_peptide_;
      Int current_mod_location_;
      ProteinHit actual_protein_;

      /// cached CV child terms for "MS:1001143"
      std::set<std::string> peptide_result_details_;

    };
  } // namespace Internal
} // namespace OpenMS
