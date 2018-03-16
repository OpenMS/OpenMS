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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MZTAB_H
#define OPENMS_FORMAT_MZTAB_H

#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <map>
#include <vector>
#include <list>
#include <algorithm>
#include <OpenMS/KERNEL/StandardTypes.h>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnon-virtual-dtor"

namespace OpenMS
{
/**
      @brief Data model of MzTab files.
      Please see the official MzTab specification at https://code.google.com/p/mztab/

      @ingroup FileIO
  */

// MzTab supports null, NaN, Inf for cells with Integer or Double values. MzTabCellType explicitly defines the state of the cell for these types.
  enum MzTabCellStateType
  {
    MZTAB_CELLSTATE_DEFAULT,
    MZTAB_CELLSTATE_NULL,
    MZTAB_CELLSTATE_NAN,
    MZTAB_CELLSTATE_INF,
    SIZE_OF_MZTAB_CELLTYPE
  };

// basic interface for all MzTab datatypes (can be null; are converted from and to cell string)
  class OPENMS_DLLAPI MzTabNullAbleInterface
  {
public:
    virtual ~MzTabNullAbleInterface();
    virtual bool isNull() const = 0;
    virtual void setNull(bool b) = 0;
    virtual String toCellString() const = 0;
    virtual void fromCellString(const String&) = 0;
  };

// interface for NaN- and Inf- able datatypes (Double and Integer in MzTab). These are as well null-able
  class OPENMS_DLLAPI MzTabNullNaNAndInfAbleInterface :
    public MzTabNullAbleInterface
  {
public:
    ~MzTabNullNaNAndInfAbleInterface() override;
    virtual bool isNaN() const = 0;
    virtual void setNaN() = 0;
    virtual bool isInf() const = 0;
    virtual void setInf() = 0;
  };

// base class for atomic, non-container types (Double, Int)
  class OPENMS_DLLAPI MzTabNullAbleBase :
    public MzTabNullAbleInterface
  {
public:
    MzTabNullAbleBase();

    ~MzTabNullAbleBase() override;

    bool isNull() const override;

    void setNull(bool b) override;

protected:
    bool null_;
  };

// base class for the atomic non-container like MzTab data types (Double, Int)
  class OPENMS_DLLAPI MzTabNullNaNAndInfAbleBase :
    public MzTabNullNaNAndInfAbleInterface
  {
public:
    MzTabNullNaNAndInfAbleBase();

    ~MzTabNullNaNAndInfAbleBase() override;

    bool isNull() const override;

    void setNull(bool b) override;

    bool isNaN() const override;

    void setNaN() override;

    bool isInf() const override;

    void setInf() override;

protected:
    MzTabCellStateType state_;
  };

  class OPENMS_DLLAPI MzTabDouble :
    public MzTabNullNaNAndInfAbleBase
  {
public:
    MzTabDouble();

    explicit MzTabDouble(const double v);

    ~MzTabDouble() override;

    void set(const double& value);

    double get() const;

    String toCellString() const override;

    void fromCellString(const String& s) override;

protected:
    double value_;
  };

  class OPENMS_DLLAPI MzTabDoubleList :
    public MzTabNullAbleBase
  {
public:
    MzTabDoubleList();

    ~MzTabDoubleList() override;

    bool isNull() const override;

    void setNull(bool b) override;

    String toCellString() const override;

    void fromCellString(const String& s) override;

    std::vector<MzTabDouble> get() const;

    void set(const std::vector<MzTabDouble>& entries);

protected:
    std::vector<MzTabDouble> entries_;
  };

  class OPENMS_DLLAPI MzTabInteger :
    public MzTabNullNaNAndInfAbleBase
  {
public:
    MzTabInteger();

    explicit MzTabInteger(const int v);

    ~MzTabInteger() override;

    void set(const Int& value);

    Int get() const;

    String toCellString() const override;

    void fromCellString(const String& s) override;

protected:
    Int value_;
  };

  class OPENMS_DLLAPI MzTabIntegerList :
    public MzTabNullAbleBase
  {
public:
    MzTabIntegerList();

    bool isNull() const override;

    void setNull(bool b) override;

    String toCellString() const override;

    void fromCellString(const String& s) override;

    std::vector<MzTabInteger> get() const;

    void set(const std::vector<MzTabInteger>& entries);

protected:
    std::vector<MzTabInteger> entries_;
  };

  class OPENMS_DLLAPI MzTabBoolean :
    public MzTabNullAbleBase
  {
public:
    MzTabBoolean();

    explicit MzTabBoolean(bool v);

    ~MzTabBoolean() override;

    void set(const bool& value);

    Int get() const;

    String toCellString() const override;

    void fromCellString(const String& s) override;

protected:
    bool value_;
  };

  class OPENMS_DLLAPI MzTabString :
    public MzTabNullAbleInterface
  {
public:
    MzTabString();

    explicit MzTabString(const String& s);

    ~MzTabString() override;

    void set(const String& value);

    String get() const;

    bool isNull() const override;

    void setNull(bool b) override;

    String toCellString() const override;

    void fromCellString(const String& s) override;

protected:
    String value_;
  };

  class OPENMS_DLLAPI MzTabParameter :
    public MzTabNullAbleInterface
  {
public:
    MzTabParameter();

    ~MzTabParameter() override;

    bool isNull() const override;

    void setNull(bool b) override;

    void setCVLabel(const String& CV_label);

    void setAccession(const String& accession);

    void setName(const String& name);

    void setValue(const String& value);

    String getCVLabel() const;

    String getAccession() const;

    String getName() const;

    String getValue() const;

    String toCellString() const override;

    void fromCellString(const String& s) override;

protected:
    String CV_label_;
    String accession_;
    String name_;
    String value_;
  };

  class OPENMS_DLLAPI MzTabParameterList :
    public MzTabNullAbleInterface
  {
public:

    ~MzTabParameterList() override;

    bool isNull() const override;

    void setNull(bool b) override;

    String toCellString() const override;

    void fromCellString(const String& s) override;

    std::vector<MzTabParameter> get() const;

    void set(const std::vector<MzTabParameter>& parameters);

protected:
    std::vector<MzTabParameter> parameters_;
  };

  class OPENMS_DLLAPI MzTabStringList :
    public MzTabNullAbleInterface
  {
public:
    MzTabStringList();

    ~MzTabStringList() override;

    // needed for e.g. ambiguity_members and GO accessions as these use ',' as separator while the others use '|'
    void setSeparator(char sep);

    bool isNull() const override;

    void setNull(bool b) override;

    String toCellString() const override;

    void fromCellString(const String& s) override;

    std::vector<MzTabString> get() const;

    void set(const std::vector<MzTabString>& entries);

protected:
    std::vector<MzTabString> entries_;
    char sep_;
  };

  class OPENMS_DLLAPI MzTabModification :
    public MzTabNullAbleInterface
  {
public:
    MzTabModification();

    ~MzTabModification() override;

    bool isNull() const override;

    void setNull(bool b) override;

    // set (potentially ambiguous) position(s) with associated parameter (might be null if not set)
    void setPositionsAndParameters(const std::vector<std::pair<Size, MzTabParameter> >& ppp);

    std::vector<std::pair<Size, MzTabParameter> > getPositionsAndParameters() const;

    void setModificationIdentifier(const MzTabString& mod_id);

    MzTabString getModOrSubstIdentifier() const;

    String toCellString() const override;

    void fromCellString(const String& s) override;

protected:
    std::vector<std::pair<Size, MzTabParameter> > pos_param_pairs_;
    MzTabString mod_identifier_;
  };

  class OPENMS_DLLAPI MzTabModificationList :
    public MzTabNullAbleBase
  {
public:
    ~MzTabModificationList() override;

    bool isNull() const override;

    void setNull(bool b) override;

    String toCellString() const override;

    void fromCellString(const String& s) override;

    std::vector<MzTabModification> get() const;

    void set(const std::vector<MzTabModification>& entries);

protected:
    std::vector<MzTabModification> entries_;

  };

  class OPENMS_DLLAPI MzTabSpectraRef :
    public MzTabNullAbleInterface
  {
public:
    MzTabSpectraRef();

    ~MzTabSpectraRef() override;

    bool isNull() const override;

    void setNull(bool b) override;

    void setMSFile(Size index);

    void setSpecRef(String spec_ref);

    String getSpecRef() const;

    Size getMSFile() const;

    void setSpecRefFile(const String& spec_ref);

    String toCellString() const override;

    void fromCellString(const String& s) override;

protected:
    Size ms_run_; // number is specified in the meta data section.
    String spec_ref_;
  };

// MTD

  struct OPENMS_DLLAPI MzTabSampleMetaData
  {
    MzTabString description;
    std::map<Size, MzTabParameter> species;
    std::map<Size, MzTabParameter> tissue;
    std::map<Size, MzTabParameter> cell_type;
    std::map<Size, MzTabParameter> disease;
    std::map<Size, MzTabParameter> custom;
  };

  struct OPENMS_DLLAPI MzTabSoftwareMetaData
  {
    MzTabParameter software;
    std::map<Size, MzTabString> setting;
  };

  struct OPENMS_DLLAPI MzTabModificationMetaData
  {
    MzTabParameter modification;
    MzTabString site;
    MzTabString position;
  };

  struct OPENMS_DLLAPI MzTabAssayMetaData
  {
    MzTabParameter quantification_reagent;
    std::map<Size, MzTabModificationMetaData> quantification_mod;
    MzTabString sample_ref;
    MzTabString ms_run_ref;
  };

  struct OPENMS_DLLAPI MzTabCVMetaData
  {
    MzTabString label;
    MzTabString full_name;
    MzTabString version;
    MzTabString url;
  };

  struct OPENMS_DLLAPI MzTabInstrumentMetaData
  {
    MzTabParameter name;
    MzTabParameter source;
    std::map<Size, MzTabParameter> analyzer;
    MzTabParameter detector;
  };

  struct OPENMS_DLLAPI MzTabContactMetaData
  {
    MzTabString name;
    MzTabString affiliation;
    MzTabString email;
  };

  struct OPENMS_DLLAPI MzTabMSRunMetaData
  {
    MzTabParameter format;
    MzTabString location;
    MzTabParameter id_format;
    MzTabParameterList fragmentation_method;
  };

  struct OPENMS_DLLAPI MzTabStudyVariableMetaData
  {
    MzTabIntegerList assay_refs;
    MzTabIntegerList sample_refs;
    MzTabString description;
  };

// all meta data of a mzTab file. Please refer to specification for documentation.
  class OPENMS_DLLAPI MzTabMetaData
  {
public:
    MzTabMetaData();

    MzTabString mz_tab_version;
    MzTabString mz_tab_mode;
    MzTabString mz_tab_type;
    MzTabString mz_tab_id;
    MzTabString title;
    MzTabString description;

    std::map<Size, MzTabParameter> protein_search_engine_score;
    std::map<Size, MzTabParameter> peptide_search_engine_score;
    std::map<Size, MzTabParameter> psm_search_engine_score;
    std::map<Size, MzTabParameter> smallmolecule_search_engine_score;

    std::map<Size, MzTabParameterList> sample_processing;

    std::map<Size, MzTabInstrumentMetaData> instrument;

    std::map<Size, MzTabSoftwareMetaData> software;

    MzTabParameterList false_discovery_rate;

    std::map<Size, MzTabString> publication;

    std::map<Size, MzTabContactMetaData> contact;

    std::map<Size, MzTabString> uri;

    std::map<Size, MzTabModificationMetaData> fixed_mod;

    std::map<Size, MzTabModificationMetaData> variable_mod;

    MzTabParameter quantification_method;

    MzTabParameter protein_quantification_unit;
    MzTabParameter peptide_quantification_unit;
    MzTabParameter small_molecule_quantification_unit;

    std::map<Size, MzTabMSRunMetaData> ms_run;

    std::map<Size, MzTabParameter> custom;

    std::map<Size, MzTabSampleMetaData> sample;

    std::map<Size, MzTabAssayMetaData> assay;

    std::map<Size, MzTabStudyVariableMetaData> study_variable;

    std::map<Size, MzTabCVMetaData> cv;

    std::vector<String> colunit_protein;
    std::vector<String> colunit_peptide;
    std::vector<String> colunit_psm;
    std::vector<String> colunit_small_molecule;
  };

  typedef std::pair<String, MzTabString> MzTabOptionalColumnEntry; //  column name (not null able), value (null able)

// PRT - Protein section (Table based)
  struct OPENMS_DLLAPI MzTabProteinSectionRow
  {
    MzTabProteinSectionRow();
    MzTabString accession; // The protein’s accession.
    MzTabString description; // Human readable description (i.e. the name)
    MzTabInteger taxid; // NEWT taxonomy for the species.
    MzTabString species; // Human readable name of the species
    MzTabString database; // Name of the protein database.
    MzTabString database_version; // String Version of the protein database.
    MzTabParameterList search_engine; // Search engine(s) identifying the protein.
    std::map<Size, MzTabDouble>  best_search_engine_score; // best_search_engine_score[1-n]
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run; // search_engine_score[index1]_ms_run[index2]
    MzTabInteger reliability;
    std::map<Size, MzTabInteger> num_psms_ms_run;
    std::map<Size, MzTabInteger> num_peptides_distinct_ms_run;
    std::map<Size, MzTabInteger> num_peptides_unique_ms_run;
    MzTabStringList ambiguity_members; // Alternative protein identifications.
    MzTabModificationList modifications; // Modifications identified in the protein.
    MzTabString uri; // Location of the protein’s source entry.
    MzTabStringList go_terms; // List of GO terms for the protein.
    MzTabDouble protein_coverage; // (0-1) Amount of protein sequence identified.
    std::map<Size, MzTabDouble> protein_abundance_assay;
    std::map<Size, MzTabDouble> protein_abundance_study_variable;
    std::map<Size, MzTabDouble> protein_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> protein_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; // Optional Columns must start with “opt_”
  };

// PEP - Peptide section (Table based)
  struct OPENMS_DLLAPI MzTabPeptideSectionRow
  {
    MzTabString sequence; // The peptide’s sequence.
    MzTabString accession; // The protein’s accession.
    MzTabBoolean unique; // 0=false, 1=true, null else: Peptide is unique for the protein.
    MzTabString database; // Name of the sequence database.
    MzTabString database_version; // Version (and optionally # of entries).
    MzTabParameterList search_engine; // Search engine(s) that identified the peptide.
    std::map<Size, MzTabDouble> best_search_engine_score; // Search engine(s) score(s) for the peptide.
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run;
    MzTabInteger reliability; // (1-3) 0=null Identification reliability for the peptide.
    MzTabModificationList modifications; // Modifications identified in the peptide.
    MzTabDoubleList retention_time; // Time points in seconds. Semantics may vary.
    MzTabDoubleList retention_time_window;
    MzTabInteger charge; // Precursor ion’s charge.
    MzTabDouble mass_to_charge; // Precursor ion’s m/z.
    MzTabString uri; // Location of the PSMs source entry.
    MzTabSpectraRef spectra_ref; // Spectra identifying the peptide.
    std::map<Size, MzTabDouble> peptide_abundance_assay;
    std::map<Size, MzTabDouble> peptide_abundance_study_variable;
    std::map<Size, MzTabDouble> peptide_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> peptide_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; // Optional columns must start with “opt_”.
  };

// PSM - PSM section (Table based)
  struct OPENMS_DLLAPI MzTabPSMSectionRow
  {
    MzTabString sequence; // The peptide’s sequence.
    MzTabInteger PSM_ID;
    MzTabString accession; // The protein’s accession.
    MzTabBoolean unique; // 0=false, 1=true, null else: Peptide is unique for the protein.
    MzTabString database; // Name of the sequence database.
    MzTabString database_version; // Version (and optionally # of entries).
    MzTabParameterList search_engine; // Search engine(s) that identified the peptide.
    std::map<Size, MzTabDouble> search_engine_score; // Search engine(s) score(s) for the peptide.
    MzTabInteger reliability; // (1-3) 0=null Identification reliability for the peptide.
    MzTabModificationList modifications; // Modifications identified in the peptide.
    MzTabDoubleList retention_time; // Time points in seconds. Semantics may vary.
    MzTabInteger charge; // The charge of the experimental precursor ion.
    MzTabDouble exp_mass_to_charge; // The m/z ratio of the experimental precursor ion.
    MzTabDouble calc_mass_to_charge;
    MzTabString uri; // Location of the PSM’s source entry.
    MzTabSpectraRef spectra_ref; // Spectra identifying the peptide.
    MzTabString pre;
    MzTabString post;
    MzTabString start;
    MzTabString end;
    std::vector<MzTabOptionalColumnEntry> opt_; // Optional columns must start with “opt_”.
  };

// SML Small molecule section (table based)
  struct OPENMS_DLLAPI MzTabSmallMoleculeSectionRow
  {
    MzTabStringList identifier; // The small molecule’s identifier.
    MzTabString chemical_formula; // Chemical formula of the identified compound.
    MzTabString smiles; // Molecular structure in SMILES format.
    MzTabString inchi_key; // InChi Key of the identified compound.
    MzTabString description; // Human readable description (i.e. the name)
    MzTabDouble exp_mass_to_charge; // Precursor ion’s m/z.
    MzTabDouble calc_mass_to_charge; // Precursor ion’s m/z.
    MzTabDouble charge; // Precursor ion’s charge.
    MzTabDoubleList retention_time; // Time points in seconds. Semantics may vary.
    MzTabInteger taxid; // NEWT taxonomy for the species.
    MzTabString species; // Human readable name of the species
    MzTabString database; // Name of the used database.
    MzTabString database_version; // String Version of the database (and optionally # of compounds).
    MzTabInteger reliability; // (1-3) The identification reliability.
    MzTabString uri; // The source entry’s location.
    MzTabSpectraRef spectra_ref; // Spectra identifying the small molecule.
    MzTabParameterList search_engine; // Search engine(s) identifying the small molecule.
    std::map<Size, MzTabDouble> best_search_engine_score; // Search engine(s) identifications score(s).
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run;
    MzTabString modifications; // Modifications identified on the small molecule.
    std::map<Size, MzTabDouble> smallmolecule_abundance_assay;
    std::map<Size, MzTabDouble> smallmolecule_abundance_study_variable;
    std::map<Size, MzTabDouble> smallmolecule_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> smallmolecule_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; // Optional columns must start with “opt_”.
  };

  typedef std::vector<MzTabProteinSectionRow> MzTabProteinSectionRows;
  typedef std::vector<MzTabPeptideSectionRow> MzTabPeptideSectionRows;
  typedef std::vector<MzTabPSMSectionRow> MzTabPSMSectionRows;
  typedef std::vector<MzTabSmallMoleculeSectionRow> MzTabSmallMoleculeSectionRows;

/**
      @brief Data model of MzTab files.
      Please see the official MzTab specification at https://code.google.com/p/mztab/

      @ingroup FileIO
 */
  class OPENMS_DLLAPI MzTab
  {
public:
    /// Default constructor
    MzTab();

    /// Destructor
    virtual ~MzTab();

    const MzTabMetaData& getMetaData() const;

    void setMetaData(const MzTabMetaData& md);

    const MzTabProteinSectionRows& getProteinSectionRows() const;

    void setProteinSectionRows(const MzTabProteinSectionRows& psd);

    const MzTabPeptideSectionRows& getPeptideSectionRows() const;

    void setPeptideSectionRows(const MzTabPeptideSectionRows& psd);

    const MzTabPSMSectionRows& getPSMSectionRows() const;

    void setPSMSectionRows(const MzTabPSMSectionRows& psd);

    void setCommentRows(const std::map<Size, String>& com);

    void setEmptyRows(const std::vector<Size>& empty);

    const std::vector<Size>& getEmptyRows() const;

    const std::map<Size, String>& getCommentRows() const;

    const MzTabSmallMoleculeSectionRows& getSmallMoleculeSectionRows() const;

    void setSmallMoleculeSectionRows(const MzTabSmallMoleculeSectionRows& smsd);

    // Extract opt_ (custom, optional column names)
    std::vector<String> getProteinOptionalColumnNames() const;

    // Extract opt_ (custom, optional column names)
    std::vector<String> getPeptideOptionalColumnNames() const;

    // Extract opt_ (custom, optional column names)
    std::vector<String> getPSMOptionalColumnNames() const;

    // Extract opt_ (custom, optional column names)
    std::vector<String> getSmallMoleculeOptionalColumnNames() const;

protected:
    MzTabMetaData meta_data_;
    MzTabProteinSectionRows protein_data_;
    MzTabPeptideSectionRows peptide_data_;
    MzTabPSMSectionRows psm_data_;
    MzTabSmallMoleculeSectionRows small_molecule_data_;
    std::vector<Size> empty_rows_; // index of empty rows
    std::map<Size, String> comment_rows_; // comments
  };

} // namespace OpenMS

#pragma clang diagnostic pop

#endif // OPENMS_FORMAT_MZTAB_H
