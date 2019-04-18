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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideEvidence.h>

#include <map>
#include <vector>
#include <list>
#include <algorithm>

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

  /// base class for the atomic non-container like MzTab data types (Double, Int)
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

    /// needed for e.g. ambiguity_members and GO accessions as these use ',' as separator while the others use '|'
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

    /// set (potentially ambiguous) position(s) with associated parameter (might be null if not set)
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
    //TODO shouldn't settings always consist of the name of the setting
    // and the value?
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
    std::vector<int> ms_run_ref; // adapted to address https://github.com/HUPO-PSI/mzTab/issues/26
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
    std::vector<int> assay_refs;
    std::vector<int> sample_refs;
    MzTabString description;
  };

  /// all meta data of a mzTab file. Please refer to specification for documentation.
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
    std::map<Size, MzTabParameter> nucleic_acid_search_engine_score;
    std::map<Size, MzTabParameter> oligonucleotide_search_engine_score;
    std::map<Size, MzTabParameter> osm_search_engine_score;

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

  typedef std::pair<String, MzTabString> MzTabOptionalColumnEntry; ///<  column name (not null able), value (null able)

  /// PRT - Protein section (Table based)
  struct OPENMS_DLLAPI MzTabProteinSectionRow
  {
    MzTabProteinSectionRow();
    MzTabString accession; ///< The protein’s accession.
    MzTabString description; ///< Human readable description (i.e. the name)
    MzTabInteger taxid; ///< NEWT taxonomy for the species.
    MzTabString species; ///< Human readable name of the species
    MzTabString database; ///< Name of the protein database.
    MzTabString database_version; ///< String Version of the protein database.
    MzTabParameterList search_engine; ///< Search engine(s) identifying the protein.
    std::map<Size, MzTabDouble>  best_search_engine_score; ///< best_search_engine_score[1-n]
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run; ///< search_engine_score[index1]_ms_run[index2]
    MzTabInteger reliability;
    std::map<Size, MzTabInteger> num_psms_ms_run;
    std::map<Size, MzTabInteger> num_peptides_distinct_ms_run;
    std::map<Size, MzTabInteger> num_peptides_unique_ms_run;
    MzTabStringList ambiguity_members; ///< Alternative protein identifications.
    MzTabModificationList modifications; ///< Modifications identified in the protein.
    MzTabString uri; ///< Location of the protein’s source entry.
    MzTabStringList go_terms; ///< List of GO terms for the protein.
    MzTabDouble coverage; ///< (0-1) Amount of protein sequence identified.
    std::map<Size, MzTabDouble> protein_abundance_assay;
    std::map<Size, MzTabDouble> protein_abundance_study_variable;
    std::map<Size, MzTabDouble> protein_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> protein_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional Columns must start with “opt_”

    /// Comparison operator for sorting rows
    struct RowCompare
    {
      bool operator()(const MzTabProteinSectionRow& row1,
                      const MzTabProteinSectionRow& row2) const
      {
        return row1.accession.get() < row2.accession.get();
      }
    };
  };

  /// PEP - Peptide section (Table based)
  struct OPENMS_DLLAPI MzTabPeptideSectionRow
  {
    MzTabString sequence; ///< The peptide’s sequence.
    MzTabString accession; ///< The protein’s accession.
    MzTabBoolean unique; ///< 0=false, 1=true, null else: Peptide is unique for the protein.
    MzTabString database; ///< Name of the sequence database.
    MzTabString database_version; ///< Version (and optionally # of entries).
    MzTabParameterList search_engine; ///< Search engine(s) that identified the peptide.
    std::map<Size, MzTabDouble> best_search_engine_score; ///< Search engine(s) score(s) for the peptide.
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run;
    MzTabInteger reliability; ///< (1-3) 0=null Identification reliability for the peptide.
    MzTabModificationList modifications; ///< Modifications identified in the peptide.
    MzTabDoubleList retention_time; ///< Time points in seconds. Semantics may vary.
    MzTabDoubleList retention_time_window;
    MzTabInteger charge; ///< Precursor ion’s charge.
    MzTabDouble mass_to_charge; ///< Precursor ion’s m/z.
    MzTabString uri; ///< Location of the PSMs source entry.
    MzTabSpectraRef spectra_ref; ///< Spectra identifying the peptide.
    std::map<Size, MzTabDouble> peptide_abundance_assay;
    std::map<Size, MzTabDouble> peptide_abundance_study_variable;
    std::map<Size, MzTabDouble> peptide_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> peptide_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.

    /// Comparison operator for sorting rows
    struct RowCompare
    {
      bool operator()(const MzTabPeptideSectionRow& row1,
                      const MzTabPeptideSectionRow& row2) const
      {
        return (std::make_pair(row1.sequence.get(), row1.accession.get()) <
                std::make_pair(row2.sequence.get(), row2.accession.get()));
      }
    };
  };

  /// PSM - PSM section (Table based)
  struct OPENMS_DLLAPI MzTabPSMSectionRow
  {
    MzTabString sequence; ///< The peptide’s sequence.
    MzTabInteger PSM_ID;
    MzTabString accession; ///< The protein’s accession.
    MzTabBoolean unique; ///< 0=false, 1=true, null else: Peptide is unique for the protein.
    MzTabString database; ///< Name of the sequence database.
    MzTabString database_version; ///< Version (and optionally # of entries).
    MzTabParameterList search_engine; ///< Search engine(s) that identified the peptide.
    std::map<Size, MzTabDouble> search_engine_score; ///< Search engine(s) score(s) for the peptide.
    MzTabInteger reliability; ///< (1-3) 0=null Identification reliability for the peptide.
    MzTabModificationList modifications; ///< Modifications identified in the peptide.
    MzTabDoubleList retention_time; ///< Time points in seconds. Semantics may vary.
    MzTabInteger charge; ///< The charge of the experimental precursor ion.
    MzTabDouble exp_mass_to_charge; ///< The m/z ratio of the experimental precursor ion.
    MzTabDouble calc_mass_to_charge;
    MzTabString uri; ///< Location of the PSM’s source entry.
    MzTabSpectraRef spectra_ref; ///< Spectra identifying the peptide.
    MzTabString pre;
    MzTabString post;
    MzTabString start;
    MzTabString end;
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.

    /// Comparison operator for sorting rows
    struct RowCompare
    {
      bool operator()(const MzTabPSMSectionRow& row1,
                      const MzTabPSMSectionRow& row2) const
      {
        // @TODO: sort by "PSM_ID"? what's the point of that field?
        return (std::make_tuple(row1.sequence.get(),
                                row1.spectra_ref.getMSFile(),
                                row1.spectra_ref.getSpecRef(),
                                row1.accession.get()) <
                std::make_tuple(row2.sequence.get(),
                                row2.spectra_ref.getMSFile(),
                                row2.spectra_ref.getSpecRef(),
                                row2.accession.get()));
      }
    };
  };

  /// SML Small molecule section (table based)
  struct OPENMS_DLLAPI MzTabSmallMoleculeSectionRow
  {
    MzTabStringList identifier; ///< The small molecule’s identifier.
    MzTabString chemical_formula; ///< Chemical formula of the identified compound.
    MzTabString smiles; ///< Molecular structure in SMILES format.
    MzTabString inchi_key; ///< InChi Key of the identified compound.
    MzTabString description; ///< Human readable description (i.e. the name)
    MzTabDouble exp_mass_to_charge; ///< Precursor ion’s m/z.
    MzTabDouble calc_mass_to_charge; ///< Precursor ion’s m/z.
    MzTabInteger charge; ///< Precursor ion’s charge.
    MzTabDoubleList retention_time; ///< Time points in seconds. Semantics may vary.
    MzTabInteger taxid; ///< NEWT taxonomy for the species.
    MzTabString species; ///< Human readable name of the species
    MzTabString database; ///< Name of the used database.
    MzTabString database_version; ///< String Version of the database (and optionally # of compounds).
    MzTabInteger reliability; ///< (1-3) The identification reliability.
    MzTabString uri; ///< The source entry’s location.
    MzTabSpectraRef spectra_ref; ///< Spectra identifying the small molecule.
    MzTabParameterList search_engine; ///< Search engine(s) identifying the small molecule.
    std::map<Size, MzTabDouble> best_search_engine_score; ///< Search engine(s) identifications score(s).
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run;
    MzTabString modifications; ///< Modifications identified on the small molecule.
    std::map<Size, MzTabDouble> smallmolecule_abundance_assay;
    std::map<Size, MzTabDouble> smallmolecule_abundance_study_variable;
    std::map<Size, MzTabDouble> smallmolecule_abundance_stdev_study_variable;
    std::map<Size, MzTabDouble> smallmolecule_abundance_std_error_study_variable;
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.
  };

  /// NUC - Nucleic acid section (table-based)
  struct OPENMS_DLLAPI MzTabNucleicAcidSectionRow
  {
    MzTabString accession; ///< The nucleic acid’s accession.
    MzTabString description; ///< Human readable description (i.e. the name)
    MzTabInteger taxid; ///< NEWT taxonomy for the species.
    MzTabString species; ///< Human readable name of the species
    MzTabString database; ///< Name of the sequence database.
    MzTabString database_version; ///< Version of the sequence database.
    MzTabParameterList search_engine; ///< Search engine(s) that identified the nucleic acid.
    std::map<Size, MzTabDouble>  best_search_engine_score; ///< Best search engine(s) score(s) (over all MS runs)
    std::map<Size, std::map<Size, MzTabDouble> > search_engine_score_ms_run;
    MzTabInteger reliability;
    std::map<Size, MzTabInteger> num_osms_ms_run;
    std::map<Size, MzTabInteger> num_oligos_distinct_ms_run;
    std::map<Size, MzTabInteger> num_oligos_unique_ms_run;
    MzTabStringList ambiguity_members; ///< Alternative nucleic acid identifications.
    MzTabModificationList modifications; ///< Modifications identified in the nucleic acid.
    MzTabString uri; ///< Location of the nucleic acid’s source entry.
    // do GO terms make sense for nucleic acid sequences?
    MzTabStringList go_terms; ///< List of GO terms for the nucleic acid.
    MzTabDouble coverage; ///< (0-1) Fraction of nucleic acid sequence identified.
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional Columns must start with “opt_”

    /// Comparison operator for sorting rows
    struct RowCompare
    {
      bool operator()(const MzTabNucleicAcidSectionRow& row1,
                      const MzTabNucleicAcidSectionRow& row2) const
      {
        return row1.accession.get() < row2.accession.get();
      }
    };
  };

  /// OLI - Oligonucleotide section (table-based)
  struct OPENMS_DLLAPI MzTabOligonucleotideSectionRow
  {
    MzTabString sequence; ///< The oligonucleotide’s sequence.
    MzTabString accession; ///< The nucleic acid’s accession.
    MzTabBoolean unique; ///< 0=false, 1=true, null else: Oligonucleotide maps uniquely to the nucleic acid sequence.
    MzTabParameterList search_engine; ///< Search engine(s) that identified the match.
    std::map<Size, MzTabDouble> best_search_engine_score; ///< Search engine(s) score(s) for the match.
    std::map<Size, std::map<Size, MzTabDouble>> search_engine_score_ms_run; ///< Search engine(s) score(s) per individual MS run
    MzTabInteger reliability; ///< (1-3) 0=null Identification reliability for the match.
    MzTabModificationList modifications; ///< Modifications identified in the oligonucleotide.
    MzTabDoubleList retention_time; ///< Time points in seconds. Semantics may vary.
    MzTabDoubleList retention_time_window;
    MzTabString uri; ///< Location of the oligonucleotide's source entry.
    MzTabString pre;
    MzTabString post;
    MzTabString start;
    MzTabString end;
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.

    /// Comparison operator for sorting rows
    struct RowCompare
    {
      bool operator()(const MzTabOligonucleotideSectionRow& row1,
                      const MzTabOligonucleotideSectionRow& row2) const
        {
          return (std::make_tuple(row1.sequence.get(), row1.accession.get(),
                                  row1.start.get(), row1.end.get()) <
                  std::make_tuple(row2.sequence.get(), row2.accession.get(),
                                  row2.start.get(), row2.end.get()));
        }
    };

  };

  /// OSM - OSM (oligonucleotide-spectrum match) section (table-based)
  struct OPENMS_DLLAPI MzTabOSMSectionRow
  {
    MzTabString sequence; ///< The oligonucleotide’s sequence.
    MzTabParameterList search_engine; ///< Search engine(s) that identified the match.
    std::map<Size, MzTabDouble> search_engine_score; ///< Search engine(s) score(s) for the match.
    MzTabInteger reliability; ///< (1-3) 0=null Identification reliability for the match.
    MzTabModificationList modifications; ///< Modifications identified in the oligonucleotide.
    MzTabDoubleList retention_time; ///< Time points in seconds. Semantics may vary.
    MzTabInteger charge; ///< The charge of the experimental precursor ion.
    MzTabDouble exp_mass_to_charge; ///< The m/z ratio of the experimental precursor ion.
    MzTabDouble calc_mass_to_charge; ///< The theoretical m/z ratio of the oligonucleotide.
    MzTabString uri; ///< Location of the OSM’s source entry.
    MzTabSpectraRef spectra_ref; ///< Reference to the spectrum underlying the match.
    std::vector<MzTabOptionalColumnEntry> opt_; ///< Optional columns must start with “opt_”.

    /// Comparison operator for sorting rows
    struct RowCompare
    {
      bool operator()(const MzTabOSMSectionRow& row1,
                      const MzTabOSMSectionRow& row2) const
      {
        return (std::make_tuple(row1.sequence.get(),
                                row1.spectra_ref.getMSFile(),
                                row1.spectra_ref.getSpecRef()) <
                std::make_tuple(row2.sequence.get(),
                                row2.spectra_ref.getMSFile(),
                                row2.spectra_ref.getSpecRef()));
      }
    };
  };

  typedef std::vector<MzTabProteinSectionRow> MzTabProteinSectionRows;
  typedef std::vector<MzTabPeptideSectionRow> MzTabPeptideSectionRows;
  typedef std::vector<MzTabPSMSectionRow> MzTabPSMSectionRows;
  typedef std::vector<MzTabSmallMoleculeSectionRow> MzTabSmallMoleculeSectionRows;
  typedef std::vector<MzTabNucleicAcidSectionRow> MzTabNucleicAcidSectionRows;
  typedef std::vector<MzTabOligonucleotideSectionRow> MzTabOligonucleotideSectionRows;
  typedef std::vector<MzTabOSMSectionRow> MzTabOSMSectionRows;


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

    const MzTabSmallMoleculeSectionRows& getSmallMoleculeSectionRows() const;

    void setSmallMoleculeSectionRows(const MzTabSmallMoleculeSectionRows& smsd);

    const MzTabNucleicAcidSectionRows& getNucleicAcidSectionRows() const;

    void setNucleicAcidSectionRows(const MzTabNucleicAcidSectionRows& nasd);

    const MzTabOligonucleotideSectionRows& getOligonucleotideSectionRows() const;

    void setOligonucleotideSectionRows(const MzTabOligonucleotideSectionRows& onsd);

    const MzTabOSMSectionRows& getOSMSectionRows() const;

    void setOSMSectionRows(const MzTabOSMSectionRows& osd);

    void setCommentRows(const std::map<Size, String>& com);

    void setEmptyRows(const std::vector<Size>& empty);

    const std::vector<Size>& getEmptyRows() const;

    const std::map<Size, String>& getCommentRows() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getProteinOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getPeptideOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getPSMOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getSmallMoleculeOptionalColumnNames() const;


    /**
      @brief Gets peptide_evidences with data from internal structures adds their info to an MzTabPSMSectionRow (pre- or unfilled)

      @param peptide_evidences Vector of PeptideEvidence holding internal data.
      @param row Pre- or unfilled MzTabPSMSectionRow to be filled with the data.
      @param rows Vector of MzTabPSMSectionRow to add the differently updated rows to.
    */
    static void addPepEvidenceToRows(const std::vector<PeptideEvidence>& peptide_evidences, MzTabPSMSectionRow& row, MzTabPSMSectionRows& rows);

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getNucleicAcidOptionalColumnNames() const;

    /// Extract opt_ (custom, optional column names)
    std::vector<String> getOligonucleotideOptionalColumnNames() const;

    static void addMetaInfoToOptionalColumns(const std::set<String>& keys, std::vector<MzTabOptionalColumnEntry>& opt, const String id, const MetaInfoInterface meta);
    /// Extract opt_ (custom, optional column names)
    std::vector<String> getOSMOptionalColumnNames() const;

    static std::map<Size, MzTabModificationMetaData> generateMzTabStringFromModifications(const std::vector<String>& mods);

    static std::map<Size, MzTabModificationMetaData> generateMzTabStringFromVariableModifications(const std::vector<String>& mods);

    static std::map<Size, MzTabModificationMetaData> generateMzTabStringFromFixedModifications(const std::vector<String>& mods);
 
    static MzTab exportFeatureMapToMzTab(const FeatureMap& feature_map, const String& filename);

    static MzTab exportIdentificationsToMzTab(
        const std::vector<ProteinIdentification>& prot_ids,
        const std::vector<PeptideIdentification>& peptide_ids,
        const String& filename,
        bool first_run_inference_only);

    /// Generate MzTab style list of PTMs from AASequence object. 
    /// All passed fixed modifications are not reported (as suggested by the standard for the PRT and PEP section).
    /// In contrast, all modifications are reported in the PSM section (see standard document for details).
    static MzTabModificationList extractModificationListFromAASequence(const AASequence& aas, const std::vector<String>& fixed_mods = std::vector<String>());

    static MzTab exportConsensusMapToMzTab(
      const ConsensusMap & consensus_map, 
      const String & filename,
      const bool export_unidentified_features,
      const bool export_unassigned_ids,
      String title = "ConsensusMap export from OpenMS");

  protected:
    /// Helper function for "get...OptionalColumnNames" functions
    template <typename SectionRows>
    std::vector<String> getOptionalColumnNames_(const SectionRows& rows) const
    {
      // vector is used to preserve the column order
      std::vector<String> names;
      if (!rows.empty())
      {
        for (typename SectionRows::const_iterator it = rows.begin(); it != rows.end(); ++it)
        {
          for (std::vector<MzTabOptionalColumnEntry>::const_iterator it_opt = it->opt_.begin(); it_opt != it->opt_.end(); ++it_opt)
          {
            if (std::find(names.begin(), names.end(), it_opt->first) == names.end())
            {
              names.push_back(it_opt->first);
            }
          }
        }
      }
      return names;
    }


    MzTabMetaData meta_data_;
    MzTabProteinSectionRows protein_data_;
    MzTabPeptideSectionRows peptide_data_;
    MzTabPSMSectionRows psm_data_;
    MzTabSmallMoleculeSectionRows small_molecule_data_;
    MzTabNucleicAcidSectionRows nucleic_acid_data_;
    MzTabOligonucleotideSectionRows oligonucleotide_data_;
    MzTabOSMSectionRows osm_data_; ///</ oligonucleotide-spectrum matches
    std::vector<Size> empty_rows_; ///< index of empty rows
    std::map<Size, String> comment_rows_; ///< comments
  };

} // namespace OpenMS

#pragma clang diagnostic pop

