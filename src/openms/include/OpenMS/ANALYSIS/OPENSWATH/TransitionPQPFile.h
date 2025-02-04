// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>

namespace OpenMS
{

  /**
      @brief This class supports reading and writing of PQP files. 

      The PQP files are SQLite databases consisting of several tables
      representing the data contained in TraML files. For another file format that stores transitions, see also
      TransitionTSVFile.

      This class can convert TraML and PQP files into each other

      <h2> The file format has the following tables: </h2>

      Genes and proteins are described by a primary key as well as a
      human-readable gene name or protein accession key.
      <table border="0"><tr><td>
        <table>
          <tr> <th BGCOLOR="#EBEBEB" colspan=3>GENE</th> </tr>
          <tr> <td BGCOLOR="#EBEBEB">ID</td> <td>INT</td> <td> Primary Key (gene id)</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">GENE_NAME</td> <td>TEXT</td> <td> Gene name </td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">DECOY</td> <td>INT (0 or 1)</td> <td> Whether this is a decoy gene (1: decoy, 0: target) </td> </tr>
        </table>
      </td><td valign="top">
        <table>
          <tr> <th BGCOLOR="#EBEBEB" colspan=3>PEPTIDE_GENE_MAPPING</th> </tr>
          <tr> <td BGCOLOR="#EBEBEB">PEPTIDE_ID</td> <td>INT</td> <td> Foreign Key (PEPTIDE.ID)</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">GENE_ID</td> <td>INT</td> <td> Foreign Key (GENE.ID)</td> </tr>
        </table>
      </td></table>

      <table border="0"><tr><td>
      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>PROTEIN</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ID</td> <td>INT</td> <td> Primary Key (protein id)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">PROTEIN_ACCESSION</td> <td>TEXT</td> <td> Protein accession </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">DECOY</td> <td>INT (0 or 1)</td> <td> Whether this is a decoy protein (1: decoy, 0: target) </td> </tr>
      </table>
      </td><td valign="top">
      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>PEPTIDE_PROTEIN_MAPPING</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">PEPTIDE_ID</td> <td>INT</td> <td> Foreign Key (PEPTIDE.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">PROTEIN_ID</td> <td>INT</td> <td> Foreign Key (PROTEIN.ID)</td> </tr>
      </table>
      </td></table>

      Peptides are physical analytes that are present in the sample and can carry post-translational modifications (PTMs). They are described by their amino-acid sequence and the modifications they carry:
      <table border="0"><tr><td>
      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>PEPTIDE</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ID</td> <td>INT</td> <td> Primary Key (peptide id)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">UNMODIFIED_SEQUENCE</td> <td>TEXT</td> <td> Peptide sequence (unmodified) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">MODIFIED_SEQUENCE</td> <td>TEXT</td> <td> Peptide sequence (modified) <sup>1</sup></td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">DECOY</td> <td>INT (0 or 1)</td> <td> Whether this is a decoy peptide (1: decoy, 0: target) </td> </tr>
      </table>
      </td><td valign="top">
      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>PRECURSOR_PEPTIDE_MAPPING</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">PRECURSOR_ID</td> <td>INT</td> <td> Foreign Key (PRECURSOR.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">PEPTIDE_ID</td> <td>INT</td> <td> Foreign Key (PEPTIDE.ID)</td> </tr>
      </table>
      </td></table>

      Compounds are generic analytes that are present in the sample (but are not peptides). This is used for small molecules which are described by their molecular formula and the SMILES representation (structural representation):
      <table border="0"><tr><td>
      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>COMPOUND</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ID</td> <td>INT</td> <td> Primary Key (compound id)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">COMPOUND_NAME</td> <td>TEXT</td> <td> Compound name (common name of the analyte)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">SUM_FORMULA</td> <td>TEXT</td> <td> Molecular formula of the compound </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">SMILES</td> <td>TEXT</td> <td> SMILES representation</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ADDUCTS</td> <td>TEXT</td> <td> List of adducts for the compound</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">DECOY</td> <td>INT (0 or 1)</td> <td> Whether this is a decoy compound (1: decoy, 0: target) </td> </tr>
      </table>
      </td><td valign="top">
      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>PRECURSOR_COMPOUND_MAPPING</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">PRECURSOR_ID</td> <td>INT</td> <td> Foreign Key (PRECURSOR.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">COMPOUND_ID</td> <td>INT</td> <td> Foreign Key (COMPOUND.ID)</td> </tr>
      </table>
      </td></table>

      Precursors are generated upon ionization from peptides or small molecule analytes (compounds) and are described by their charge state, mass-to-charge ratio, retention time and ion mobility drift time:
      <table border="0"><tr><td>
      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>PRECURSOR</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ID</td> <td>INT</td> <td> Primary Key (precursor id)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TRAML_ID</td> <td>TEXT</td> <td> TraML identifiers (maps to the "id=" attribute in TraML)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">GROUP_LABEL</td> <td>TEXT</td> <td> designates to which peptide label group (as defined in MS:1000893) the peptide belongs to<sup>2</sup></td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">PRECURSOR_MZ</td> <td>REAL</td> <td> %Precursor m/z </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">CHARGE</td> <td>TEXT</td> <td> %Precursor charge state </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">LIBRARY_INTENSITY</td> <td>TEXT</td> <td> %Precursor library intensity </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">LIBRARY_RT</td> <td>TEXT</td> <td> Library retention time (RT) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">LIBRARY_DRIFT_TIME</td> <td>TEXT</td> <td> Library drift time (ion mobility drift time or collisional cross-section) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">DECOY</td> <td>INT (0 or 1)</td> <td> Whether this is a decoy precursor (1: decoy, 0: target) </td> </tr>
      </table>
      </td><td valign="top">
      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>TRANSITION_PRECURSOR_MAPPING</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TRANSITION_ID</td> <td>INT</td> <td> Foreign Key (TRANSITION.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">PRECURSOR_ID</td> <td>INT</td> <td> Foreign Key (PRECURSOR.ID)</td> </tr>
      </table>
      </td></table>


      Transitions are generated upon fragmentation from precursors and are described by their charge state and (fragment) mass-to-charge ratio. For peptide fragments, an ion type (e.g. y for y-ions) and a ordinal (e.g. 6 for a y6 ion) can be recorded. Note that detecting transitions are used for precursor detection and will indicate whether a given precursor is present or not in the sample. Use detecting transitions for the top N transitions of a precursor to detect a set of transitions in the sample. Identifying transitions will be used to discriminate different peptidoforms of the same precursor (see <a href="http://openswath.org/en/latest/docs/ipf.html">IPF Workflow</a> for PTM inference).
      <table border="0"><tr><td>
      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>TRANSITION</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ID</td> <td>INT</td> <td> Primary Key (transition id)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TRAML_ID</td> <td>TEXT</td> <td> TraML identifiers (maps to the "id=" attribute in TraML)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">PRODUCT_MZ</td> <td>TEXT</td> <td> Fragment ion m/z </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">CHARGE</td> <td>TEXT</td> <td> Fragment ion charge </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TYPE</td> <td>CHAR</td> <td> Fragment ion type (e.g. "b" or "y")</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ANNOTATION</td> <td>TEXT</td> <td> Fragment ion annotation </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ORDINAL</td> <td>INT</td> <td> Fragment ion ordinal </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">DETECTING</td> <td>INT (0 or 1)</td> <td>1: use transition to detect peak group, 0: don't use transition for detection</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">IDENTIFYING</td> <td>INT (0 or 1)</td> <td> 1: use transition for peptidoform inference in the <a href="http://openswath.org/en/latest/docs/ipf.html">IPF Workflow</a>, 0: don't use transition for identification</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">QUANTIFYING</td> <td>INT (0 or 1)</td> <td> 1: use transition to quantify peak group, 0: don't use transition for quantification</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">LIBRARY_INTENSITY</td> <td>REAL</td> <td> Fragment ion library intensity </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">DECOY</td> <td>INT (0 or 1)</td> <td> Whether this is a decoy transition (1: decoy, 0: target) </td> </tr>
      </table>
      </td>
      </table>
      </td></table>

      There is one extra table directly mapping TRANSITION to PEPTIDE which is mainly used for the 
      <a href="http://openswath.org/en/latest/docs/ipf.html">IPF Workflow</a> for PTM inference. It directly maps transitions to peptidoforms 
      (one identification transition can map to multiple peptidoforms):

      <table border="0"><tr><td>
      </td><td>
      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>TRANSITION_PEPTIDE_MAPPING</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TRANSITION_ID</td> <td>INT</td> <td> Foreign Key (TRANSITION.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">PEPTIDE_ID</td> <td>INT</td> <td> Foreign Key (PEPTIDE.ID)</td> </tr>
      </table>
      </td></table>
      </td></table>

      Another extra table describes the file format version:
      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>VERSION</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ID</td> <td>INT</td> <td> %File Format version </td> </tr>
      </table>

      <p>
      Remarks:
      </p>
      <ul>
        <li>
          1. modifications should be supplied inside the sequence using UniMod
            identifiers or freetext identifiers that are understood by %OpenMS. See also @ref OpenMS::AASequence for more information. For example:
            <ul>
            <li> PEPT(Phosphorylation)IDE(UniMod:27)A ) </li>
            </ul>
        </li>
        <li>
          2. peptide label groups designate groups of peptides that are isotopically
          modified forms of the same peptide species. For example, the heavy and
          light forms of the same peptide will both be assigned the same peptide
          group label. For example:
            <ul>
            <li> PEPTIDEAK -> gets label "PEPTIDEAK_gr1"  </li>
            <li> PEPTIDEAK[+8] -> gets label "PEPTIDEAK_gr1"  </li>
            <li> PEPT(Phosphorylation)IDEAK -> gets label "PEPTIDEAK_gr2"  </li>
            <li> PEPT(Phosphorylation)IDEAK[+8] -> gets label "PEPTIDEAK_gr2"  </li>
            </ul>
        </li>
      </ul>
      </p>


      @htmlinclude OpenMS_TransitionPQPFile.parameters

  */
  class OPENMS_DLLAPI TransitionPQPFile :
    public TransitionTSVFile
  {

private:

    /** @brief Read PQP SQLite file
     *
     * @param filename The input file
     * @param transition_list The output list of transitions
     * @param legacy_traml_id Should legacy TraML IDs be used (boolean)?
     *
    */
    void readPQPInput_(const char* filename, std::vector<TSVTransition>& transition_list, bool legacy_traml_id = false);

    /** @brief Write a TargetedExperiment to a file
     *
     * @param filename Name of the output file
     * @param targeted_exp The data structure to be written to the file
    */
    void writePQPOutput_(const char* filename, OpenMS::TargetedExperiment& targeted_exp);

public:

    //@{
    /// Constructor
    TransitionPQPFile();

    /// Destructor
    ~TransitionPQPFile() override;
    //@}

    /** @brief Write out a targeted experiment (TraML structure) into a PQP file
     *
      @param filename The output file
      @param targeted_exp The targeted experiment
     *
    */
    void convertTargetedExperimentToPQP(const char* filename, OpenMS::TargetedExperiment& targeted_exp);

    /** @brief Read in a PQP file and construct a targeted experiment (TraML structure)
     *
      @param filename The input file
      @param targeted_exp The output targeted experiment
      @param legacy_traml_id Should legacy TraML IDs be used (boolean)?
     *
    */
    void convertPQPToTargetedExperiment(const char* filename, OpenMS::TargetedExperiment& targeted_exp, bool legacy_traml_id = false);

    /** @brief Read in a PQP file and construct a targeted experiment (Light transition structure)
     *
     * @param filename The input file
     * @param targeted_exp The output targeted experiment
     * @param legacy_traml_id Should legacy TraML IDs be used (boolean)?
     *
    */
    void convertPQPToTargetedExperiment(const char* filename, OpenSwath::LightTargetedExperiment& targeted_exp, bool legacy_traml_id = false);

    /** @brief Creates an undordered map between the traml_id and the pqp id
     * 
     * @param filename The input file
     * @param tableName The name of the table (can be "PRECURSOR" or "TRANSITION" since theses are the only tables that have a TRAML_ID)
     */
    std::unordered_map<std::string, std::string> getPQPIDToTraMLIDMap(const char* filename, std::string tableName);

  };
}