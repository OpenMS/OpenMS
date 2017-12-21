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
// $Maintainer:  $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow, Marc Sturm, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>

#include <boost/math/special_functions/fpclassify.hpp>

#include <vector>
#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_TextExporter TextExporter

  @brief This application converts several %OpenMS XML formats (featureXML, consensusXML, and idXML) to text files.

  <CENTER>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ TextExporter \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> almost any TOPP tool </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> external tools (MS Excel, OpenOffice, Notepad)</td>
  </tr>
  </table>
  </CENTER>

  The goal of this tool is to create output in a table format that is easily readable in Excel or OpenOffice. Lines in the output correspond to rows in the table; the individual columns are delineated by a separator, e.g. tab (default, TSV format) or comma (CSV format).

  Output files begin with comment lines, starting with the special character "#". The last such line(s) will be a header with column names, but this may be preceded by more general comments.

  Because the OpenMS XML formats contain different kinds of data in a hierarchical structure, TextExporter produces somewhat unusual TSV/CSV files for many inputs: Different lines in the output may belong to different types of data, and the number of columns and the meanings of the individual fields depend on the type. In such cases, the first column always contains an indicator (in capital letters) for the data type of the current line. In addition, some lines have to be understood relative to a previous line, if there is a hierarchical relationship in the data. (See below for details and examples.)

  Missing values are represented by "-1" or "nan" in numeric fields and by blanks in character/text fields.

  Depending on the input and the parameters, the output contains the following columns:

  <B>featureXML input:</B>
  - first column: @p RUN / @p PROTEIN / @p UNASSIGNEDPEPTIDE / @p FEATURE / @p PEPTIDE (indicator for the type of data in the current row)
  - a @p RUN line contains information about a protein identification run; further columns: @p run_id, @p score_type, @p score_direction, @p data_time, @p search_engine_version, @p parameters
  - a @p PROTEIN line contains data of a protein identified in the previously listed run; further columns: @p score, @p rank, @p accession, @p coverage, @p sequence
  - an @p UNASSIGNEDPEPTIDE line contains data of peptide hit that was not assigned to any feature; further columns: @p rt, @p mz, @p score, @p rank, @p sequence, @p charge, @p aa_before, @p aa_after, @p score_type, @p search_identifier, @p accessions
  - a @p FEATURE line contains data of a single feature; further columns: @p rt, @p mz, @p intensity, @p charge, @p width, @p quality, @p rt_quality, @p mz_quality, @p rt_start, @p rt_end
  - a @p PEPTIDE line contains data of a peptide hit annotated to the previous feature; further columns: same as for @p UNASSIGNEDPEPTIDE

  With the @p no_ids flag, only @p FEATURE lines (without the @p FEATURE indicator) are written.

  With the @p feature:minimal flag, only the @p rt, @p mz, and @p intensity columns of @p FEATURE lines are written.

  <B>consensusXML input:</B>

  Output format produced for the @p out parameter:
  - first column: @p MAP / @p RUN / @p PROTEIN / @p UNASSIGNEDPEPTIDE / @p CONSENSUS / @p PEPTIDE (indicator for the type of data in the current row)
  - a @p MAP line contains information about a sub-map; further columns: @p id, @p filename, @p label, @p size (potentially followed by further columns containing meta data, depending on the input)
  - a @p CONSENSUS line contains data of a single consensus feature; further columns: @p rt_cf, @p mz_cf, @p intensity_cf, @p charge_cf, @p width_cf, @p quality_cf, @p rt_X0, @p mz_X0, ..., rt_X1, mz_X1, ...
  - @p "..._cf" columns refer to the consensus feature itself, @p "..._Xi" columns refer to a sub-feature from the map with ID "Xi" (no @p quality column in this case); missing sub-features are indicated by "nan" values
  - see above for the formats of @p RUN, @p PROTEIN, @p UNASSIGNEDPEPTIDE, @p PEPTIDE lines

  With the @p no_ids flag, only @p MAP and @p CONSENSUS lines are written.

  Output format produced for the @p consensus_centroids parameter:
  - one line per consensus centroid
  - columns: @p rt, @p mz, @p intensity, @p charge, @p width, @p quality

  Output format produced for the @p consensus_elements parameter:
  - one line per sub-feature (element) of a consensus feature
  - first column: @p H / @p L (indicator for new/repeated element)
  - @p H indicates a new element, @p L indicates the replication of the first element of the current consensus feature (for plotting)
  - further columns: @p rt, @p mz, @p intensity, @p charge, @p width, @p rt_cf, @p mz_cf, @p intensity_cf, @p charge_cf, @p width_cf, @p quality_cf
  - @p "..._cf" columns refer to the consensus feature, the other columns refer to the sub-feature

  Output format produced for the @p consensus_features parameter:
  - one line per consensus feature (suitable for processing with e.g. <a href="http://www.r-project.org">R</a>)
  - columns: same as for a @p CONSENSUS line above, followed by additional columns for identification data
  - additional columns: @p peptide_N0, @p n_diff_peptides_N0, @p protein_N0, @p n_diff_proteins_N0, @p peptide_N1, ...
  - @p "..._Ni" columns refer to the identification run with index "Ni", @p n_diff_... stands for "number of different ..."; different peptides/proteins in one column are separated by "/"

  With the @p no_ids flag, the additional columns are not included.

  <B>idXML input:</B>
  - first column: @p RUN / @p PROTEIN / @p PEPTIDE (indicator for the type of data in the current row)
  - see above for the formats of @p RUN, @p PROTEIN, @p PEPTIDE lines
  - additional column for @p PEPTIDE lines: @p predicted_rt (predicted retention time)
  - additional column for @p PEPTIDE lines: @p predicted_pt (predicted proteotypicity)

  With the @p id:proteins_only flag, only @p RUN and @p PROTEIN lines are written.

  With the @p id:peptides_only flag, only @p PEPTIDE lines (without the @p PEPTIDE indicator) are written.

  With the @p id:first_dim_rt flag, the additional columns @p rt_first_dim and @p predicted_rt_first_dim are included for @p PEPTIDE lines.

  @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_TextExporter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_TextExporter.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{
  // write data from a feature to the output stream
  void writeFeature(SVOutStream& out, Peak2D::CoordinateType rt,
                    Peak2D::CoordinateType mz, Peak2D::IntensityType intensity,
                    Int charge, BaseFeature::WidthType width)
  {
    out.writeValueOrNan(rt);
    out.writeValueOrNan(mz);
    out.writeValueOrNan(intensity);
    out << charge;
    out.writeValueOrNan(width);
  }

  // stream output operator for FeatureHandle
  SVOutStream& operator<<(SVOutStream& out, const FeatureHandle& feature)
  {
    writeFeature(out, feature.getRT(), feature.getMZ(), feature.getIntensity(),
                 feature.getCharge(), feature.getWidth());
    return out;
  }

  // general stream output operator for features and consensus features
  SVOutStream& operator<<(SVOutStream& out, const BaseFeature& feature)
  {
    writeFeature(out, feature.getRT(), feature.getMZ(), feature.getIntensity(),
                 feature.getCharge(), feature.getWidth());
    out.writeValueOrNan(feature.getQuality());
    return out;
  }

  // stream output operator for consensus features
  SVOutStream& operator<<(SVOutStream& out, const ConsensusFeature& feature)
  {
    return out << static_cast<const BaseFeature&>(feature);
  }

  // stream output operator for features
  SVOutStream& operator<<(SVOutStream& out, const Feature& feature)
  {
    return out << static_cast<const BaseFeature&>(feature);
  }

  // write the header for feature data
  void writeFeatureHeader(SVOutStream& out, const String& suffix = "",
                          bool incl_quality = true, bool comment = true)
  {
    StringList elements = ListUtils::create<String>("#rt,mz,intensity,charge,width");
    if (!comment) elements[0] = "rt";
    if (incl_quality) elements.push_back("quality");
    bool old = out.modifyStrings(false);
    for (StringList::iterator it = elements.begin(); it != elements.end(); ++it)
    {
      out << *it + suffix;
    }
    out.modifyStrings(old);
  }

  // write the header for exporting consensusXML
  void writeConsensusHeader(SVOutStream& out, const String& what,
                            const String& infile, const String& now,
                            const StringList& add_comments = StringList())
  {
    out.write("#" + what + " extracted from " + infile + " on " + now + "\n");
    for (StringList::const_iterator it = add_comments.begin();
         it != add_comments.end(); ++it)
    {
      out.write("#" + *it + "\n");
    }
  }

  // write the header for run data
  void writeRunHeader(SVOutStream& out)
  {
    bool old = out.modifyStrings(false);
    out << "#RUN" << "run_id" << "score_type" << "score_direction"
        << "date_time" << "search_engine_version" << "parameters" << nl;
    out.modifyStrings(old);
  }

  // write the header for protein data
  void writeProteinHeader(SVOutStream& out)
  {
    bool old = out.modifyStrings(false);
    out << "#PROTEIN" << "score" << "rank" << "accession" << "protein_description" << "coverage"
        << "sequence" << nl;
    out.modifyStrings(old);
  }


  void writeMetaValuesHeader(SVOutStream& output, const StringList& meta_keys)
  {
    if (!meta_keys.empty())
    {
      for (StringList::const_iterator its = meta_keys.begin(); its != meta_keys.end(); ++its)
      {
        output << *its;
      }
    }
  }

  template<typename T>
  void writeMetaValues(SVOutStream& output, const T& meta_value_provider, const StringList& meta_keys)
  {
    if (!meta_keys.empty())
    {
      for (StringList::const_iterator its = meta_keys.begin(); its != meta_keys.end(); ++its)
      {
        if (meta_value_provider.metaValueExists(*its))
        {
          output << meta_value_provider.getMetaValue(*its);
        }
        else
        {
          output << "";
        }
      }
    }
  }

  // stream output operator for a ProteinHit
  SVOutStream& operator<<(SVOutStream& out, const ProteinHit& hit)
  {    
    out << hit.getScore() << hit.getRank() << hit.getAccession() << hit.getDescription()
        << hit.getCoverage() << hit.getSequence();
    return out;
  }

  // stream output operator for SearchParameters
  SVOutStream& operator<<(SVOutStream& out,
                          const ProteinIdentification::SearchParameters sp)
  {
    String param_line = "db=" + sp.db + ", db_version=" +   sp.db_version +
                        ", taxonomy=" + sp.taxonomy + ", charges=" + sp.charges + ", mass_type=";
    if (sp.mass_type == ProteinIdentification::MONOISOTOPIC)
    {
      param_line += "monoisotopic";
    }
    else param_line += "average";
    param_line += ", fixed_modifications=";
    for (vector<String>::const_iterator mit = sp.fixed_modifications.begin();
         mit != sp.fixed_modifications.end(); ++mit)
    {
      if (mit != sp.fixed_modifications.begin())
      {
        param_line += ";";
      }
      param_line += *mit;
    }
    param_line += ", variable_modifications=";
    for (vector<String>::const_iterator mit = sp.variable_modifications.begin();
         mit != sp.variable_modifications.end(); ++mit)
    {
      if (mit != sp.variable_modifications.begin())
      {
        param_line += ";";
      }
      param_line += *mit;
    }
    param_line += ", enzyme=";
    param_line += sp.digestion_enzyme.getName();
    param_line += ", missed_cleavages=" + String(sp.missed_cleavages) +
                  ", peak_mass_tolerance=" + String(sp.fragment_mass_tolerance) +
                  ", precursor_mass_tolerance=" + String(sp.precursor_mass_tolerance);
    out << param_line;
    return out;
  }

  // write a protein identification to the output stream
  void writeProteinId(SVOutStream& out, const ProteinIdentification& pid)
  {
    // protein id header
    out << "RUN" << pid.getIdentifier() << pid.getScoreType();
    if (pid.isHigherScoreBetter()) out << "higher-score-better";
    else out << "lower-score-better";
    // using ISODate ensures that TOPP tests will run through regardless of
    // locale setting
    out << pid.getDateTime().toString(Qt::ISODate).toStdString()
        << pid.getSearchEngineVersion();
    // search parameters
    ProteinIdentification::SearchParameters sp = pid.getSearchParameters();
    out << sp << nl;
    for (vector<ProteinHit>::const_iterator hit_it = pid.getHits().begin();
         hit_it != pid.getHits().end(); ++hit_it)
    {
      out << "PROTEIN" << *hit_it << nl;
    }
  }

  // write the header for peptide data
  void writePeptideHeader(SVOutStream& out, const String& what = "PEPTIDE",
                          bool incl_pred_rt = false,
                          bool incl_pred_pt = false,
                          bool incl_first_dim = false)
  {
    bool old = out.modifyStrings(false);
    if (what.empty()) out << "#rt";
    else out << "#" + what << "rt";
    out << "mz" << "score" << "rank" << "sequence" << "charge" << "aa_before"
        << "aa_after" << "score_type" << "search_identifier" << "accessions";
    if (incl_pred_rt) out << "predicted_rt";

    if (incl_first_dim) out << "rt_first_dim" << "predicted_rt_first_dim";

    if (incl_pred_pt) out << "predicted_pt";

    out.modifyStrings(old);
  }

  // stream output operator for a PeptideHit
  // TODO: output of multiple peptide evidences
  SVOutStream& operator<<(SVOutStream& out, const PeptideHit& hit)
  {
    vector<PeptideEvidence> pes = hit.getPeptideEvidences();

    if (!pes.empty())
    {
      out << hit.getScore() << hit.getRank() << hit.getSequence()
          << hit.getCharge() << pes[0].getAABefore() << pes[0].getAAAfter();
    }
    else
    {
      out << hit.getScore() << hit.getRank() << hit.getSequence()
          << hit.getCharge() << PeptideEvidence::UNKNOWN_AA << PeptideEvidence::UNKNOWN_AA;
    }
    return out;
  }

  // write a peptide identification to the output stream
  void writePeptideId(SVOutStream& out, const PeptideIdentification& pid,
                      const String& what = "PEPTIDE", bool incl_pred_rt = false, bool incl_pred_pt = false,
                      bool incl_first_dim = false, const StringList& peptide_id_meta_keys = StringList(), const StringList& peptide_hit_meta_keys = StringList())
  {
    for (vector<PeptideHit>::const_iterator hit_it = pid.getHits().begin();
         hit_it != pid.getHits().end(); ++hit_it)
    {
      if (!what.empty())
      {
        out << what;
      }

      if (pid.hasRT())
      {
        out << pid.getRT();
      }
      else
      {
        out << "-1";
      }

      if (pid.hasMZ())
      {
        out << pid.getMZ();
      }
      else
      {
        out << "-1";
      }

      out << *hit_it << pid.getScoreType() << pid.getIdentifier();

      String accessions;
      set<String> protein_accessions = hit_it->extractProteinAccessionsSet();
      for (set<String>::const_iterator acc_it = protein_accessions.begin(); acc_it != protein_accessions.end(); ++acc_it)
      {
        if (acc_it != protein_accessions.begin())
        {
          accessions += ";";
        }
        accessions += *acc_it;
      }
      out << accessions;

      if (incl_pred_rt)
      {
        if (hit_it->metaValueExists("predicted_RT"))
        {
          out << hit_it->getMetaValue("predicted_RT");
        }
        else out << "-1";
      }
      if (incl_first_dim)
      {
        if (pid.metaValueExists("first_dim_rt"))
        {
          out << pid.getMetaValue("first_dim_rt");
        }
        else out << "-1";
        if (hit_it->metaValueExists("predicted_RT_first_dim"))
        {
          out << hit_it->getMetaValue("predicted_RT_first_dim");
        }
        else out << "-1";
      }
      if (incl_pred_pt)
      {
        if (hit_it->metaValueExists("predicted_PT"))
        {
          out << hit_it->getMetaValue("predicted_PT");
        }
        else out << "-1";
      }
      writeMetaValues(out, pid, peptide_id_meta_keys);
      writeMetaValues(out, *hit_it, peptide_hit_meta_keys);
      out << nl;
    }
  }

  class TOPPTextExporter :
    public TOPPBase
  {
public:
    TOPPTextExporter() :
      TOPPBase("TextExporter", "Exports various XML formats to a text file.")
    {
    }

protected:

    void registerOptionsAndFlags_() override
    {
      registerInputFile_("in", "<file>", "", "Input file ");
      setValidFormats_("in", ListUtils::create<String>("featureXML,consensusXML,idXML,mzML"));
      registerOutputFile_("out", "<file>", "", "Output file (mandatory for featureXML and idXML)", false);
      setValidFormats_("out", ListUtils::create<String>("csv"));
      registerStringOption_("separator", "<sep>", "", "The used separator character(s); if not set the 'tab' character is used", false);
      registerStringOption_("replacement", "<string>", "_", "Used to replace occurrences of the separator in strings before writing, if 'quoting' is 'none'", false);
      registerStringOption_("quoting", "<method>", "none", "Method for quoting of strings: 'none' for no quoting, 'double' for quoting with doubling of embedded quotes,\n'escape' for quoting with backslash-escaping of embedded quotes", false);
      setValidStrings_("quoting", ListUtils::create<String>("none,double,escape"));
      registerFlag_("no_ids", "Suppresses output of identification data.");
      addEmptyLine_();

      registerTOPPSubsection_("feature", "Options for featureXML input files");
      registerFlag_("feature:minimal", "Set this flag to write only three attributes: RT, m/z, and intensity.");
      registerIntOption_("feature:add_metavalues", "<min_frequency>", -1, "Add columns for meta values which occur with a certain frequency (0-100%). Set to -1 to omit meta values (default).", false);
      setMinInt_("feature:add_metavalues", -1);
      setMaxInt_("feature:add_metavalues", 100);
      addEmptyLine_();

      registerTOPPSubsection_("id", "Options for idXML input files");
      registerFlag_("id:proteins_only", "Set this flag if you want only protein information from an idXML file");
      registerFlag_("id:peptides_only", "Set this flag if you want only peptide information from an idXML file");
      registerFlag_("id:first_dim_rt", "If this flag is set the first_dim RT of the peptide hits will also be printed (if present).");
      registerIntOption_("id:add_metavalues", "<min_frequency>", -1, "Add columns for meta values which occur with a certain frequency (0-100%). Set to -1 to omit meta values (default).", false);
      setMinInt_("id:add_metavalues", -1);
      setMaxInt_("id:add_metavalues", 100);
      registerIntOption_("id:add_hit_metavalues", "<min_frequency>", -1, "Add columns for meta values which occur with a certain frequency (0-100%). Set to -1 to omit meta values (default).", false);
      setMinInt_("id:add_hit_metavalues", -1);
      setMaxInt_("id:add_hit_metavalues", 100);
      addEmptyLine_();

      registerTOPPSubsection_("consensus", "Options for consensusXML input files");
      registerOutputFile_("consensus:centroids", "<file>", "", "Output file for centroids of consensus features", false);
      setValidFormats_("consensus:centroids", ListUtils::create<String>("csv"));
      registerOutputFile_("consensus:elements", "<file>", "", "Output file for elements of consensus features", false);
      setValidFormats_("consensus:elements", ListUtils::create<String>("csv"));
      registerOutputFile_("consensus:features", "<file>", "", "Output file for consensus features and contained elements from all maps (writes 'nan's if elements are missing)", false);
      setValidFormats_("consensus:features", ListUtils::create<String>("csv"));
      registerStringOption_("consensus:sorting_method", "<method>", "none", "Sorting options can be combined. The precedence is: sort_by_size, sort_by_maps, sorting_method", false);
      setValidStrings_("consensus:sorting_method", ListUtils::create<String>("none,RT,MZ,RT_then_MZ,intensity,quality_decreasing,quality_increasing"));
      registerFlag_("consensus:sort_by_maps", "Apply a stable sort by the covered maps, lexicographically", false);
      registerFlag_("consensus:sort_by_size", "Apply a stable sort by decreasing size (i.e., the number of elements)", false);
    }

    ExitCodes main_(int, const char**) override
    {
      //-------------------------------------------------------------
      // parameter handling
      //-------------------------------------------------------------
      String in = getStringOption_("in");
      String out = getStringOption_("out");
      bool no_ids = getFlag_("no_ids");
      bool first_dim_rt = getFlag_("id:first_dim_rt");
      int add_feature_metavalues = getIntOption_("feature:add_metavalues");
      int add_id_metavalues = getIntOption_("id:add_metavalues");
      int add_hit_metavalues = getIntOption_("id:add_hit_metavalues");

      // separator etc.
      String sep = getStringOption_("separator");
      if (sep == "") sep = "\t";
      String replacement = getStringOption_("replacement");
      String quoting = getStringOption_("quoting");
      String::QuotingMethod quoting_method;
      if (quoting == "none") quoting_method = String::NONE;
      else if (quoting == "double") quoting_method = String::DOUBLE;
      else quoting_method = String::ESCAPE;

      // input file type
      FileTypes::Type in_type = FileHandler::getType(in);
      writeDebug_(String("Input file type: ") +
                  FileTypes::typeToName(in_type), 2);

      if (in_type == FileTypes::UNKNOWN)
      {
        writeLog_("Error: Could not determine input file type!");
        return PARSE_ERROR;
      }

      StringList meta_keys;

      if (in_type == FileTypes::FEATUREXML)
      {
        //-------------------------------------------------------------
        // loading input
        //-------------------------------------------------------------

        FeatureMap feature_map;
        FeatureXMLFile f;
        f.load(in, feature_map);

        // extract common id and hit meta values
        StringList peptide_id_meta_keys;
        StringList peptide_hit_meta_keys;

        vector<PeptideIdentification> pids;
        if (add_id_metavalues >= 0 || add_hit_metavalues >= 0)
        {
                const vector<PeptideIdentification>& uapids = feature_map.getUnassignedPeptideIdentifications();
                pids.insert(pids.end(), uapids.begin(), uapids.end());
                for (FeatureMap::const_iterator cmit = feature_map.begin(); cmit != feature_map.end(); ++cmit)
                {
                        const vector<PeptideIdentification>& cpids = cmit->getPeptideIdentifications();
                        pids.insert(pids.end(), cpids.begin(), cpids.end());
                }
                if (add_id_metavalues >= 0)
                {
                        peptide_id_meta_keys = MetaInfoInterfaceUtils::findCommonMetaKeys<vector<PeptideIdentification>, StringList>(pids.begin(), pids.end(), add_id_metavalues);
                        // currently there is some hardcoded logic to create extra columns for these meta values so remove them to prevent duplication
                        peptide_id_meta_keys.erase(std::remove(peptide_id_meta_keys.begin(), peptide_id_meta_keys.end(), "predicted_RT"), peptide_id_meta_keys.end());
                        peptide_id_meta_keys.erase(std::remove(peptide_id_meta_keys.begin(), peptide_id_meta_keys.end(), "predicted_RT_first_dim"), peptide_id_meta_keys.end());
                        peptide_id_meta_keys.erase(std::remove(peptide_id_meta_keys.begin(), peptide_id_meta_keys.end(), "first_dim_rt"), peptide_id_meta_keys.end());
                        peptide_id_meta_keys.erase(std::remove(peptide_id_meta_keys.begin(), peptide_id_meta_keys.end(), "predicted_PT"), peptide_id_meta_keys.end());
                }
                if (add_hit_metavalues >= 0)
                {
                        vector<PeptideHit> temp_hits;
                        for (Size i = 0; i != pids.size(); ++i)
                        {
                                const vector<PeptideHit>& hits = pids[i].getHits();
                                temp_hits.insert(temp_hits.end(), hits.begin(), hits.end());
                        }

                        // siehe oben / analog machen
                        peptide_hit_meta_keys = MetaInfoInterfaceUtils::findCommonMetaKeys<vector<PeptideHit>, StringList>(temp_hits.begin(), temp_hits.end(), add_hit_metavalues);
                }
        }

        if (add_feature_metavalues >= 0) 
        {
          meta_keys = MetaInfoInterfaceUtils::findCommonMetaKeys<FeatureMap, StringList>(feature_map.begin(), feature_map.end(), add_feature_metavalues);
        }

        vector<ProteinIdentification> prot_ids = feature_map.getProteinIdentifications();

        // text output
        ofstream outstr(out.c_str());
        SVOutStream output(outstr, sep, replacement, quoting_method);

        bool minimal = getFlag_("feature:minimal");
        no_ids |= minimal; // "minimal" implies "no_ids"

        // write header:
        output.modifyStrings(false);
        bool comment = true;
        if (!no_ids)
        {
          writeRunHeader(output);
          writeProteinHeader(output);
          writePeptideHeader(output, "UNASSIGNEDPEPTIDE");
          writeMetaValuesHeader(output, peptide_id_meta_keys);
          writeMetaValuesHeader(output, peptide_hit_meta_keys);
          output << nl;
          output << "#FEATURE";
          comment = false;
        }
        if (minimal) output << "#rt" << "mz" << "intensity";
        else
        {
          writeFeatureHeader(output, "", true, comment);
          output << "rt_quality" << "mz_quality" << "rt_start" << "rt_end";
        }
        writeMetaValuesHeader(output, meta_keys);
        output << nl;
        if (!no_ids)
        {
          writePeptideHeader(output);
          writeMetaValuesHeader(output, peptide_id_meta_keys);
          writeMetaValuesHeader(output, peptide_hit_meta_keys);
          output << nl;
        }
        output.modifyStrings(true);

        if (!no_ids)
        {
          for (vector<ProteinIdentification>::const_iterator it =
                 prot_ids.begin(); it != prot_ids.end(); ++it)
          {
            writeProteinId(output, *it);
          }
          for (vector<PeptideIdentification>::const_iterator pit =
                 feature_map.getUnassignedPeptideIdentifications().begin();
               pit != feature_map.getUnassignedPeptideIdentifications().end();
               ++pit)
          {
            writePeptideId(output, *pit, "UNASSIGNEDPEPTIDE", false, false, false, peptide_id_meta_keys, peptide_hit_meta_keys);
          }
        }

        for (FeatureMap::const_iterator citer = feature_map.begin();
             citer != feature_map.end(); ++citer)
        {
          if (!no_ids)
          {
            output << "FEATURE";
          }
          if (minimal)
          {
            output << citer->getRT() << citer->getMZ()
                   << citer->getIntensity();
          }
          else
          {
            output << *citer << citer->getQuality(0) << citer->getQuality(1);
            if (citer->getConvexHulls().size() > 0)
            {
              output << citer->getConvexHulls().begin()->getBoundingBox().minX() 
                     << citer->getConvexHulls().begin()->getBoundingBox().maxX();
            }
            else
            {
              output << "-1" << "-1";
            }
          }
          writeMetaValues(output, *citer, meta_keys);
          output << nl;

          // peptide ids
          if (!no_ids)
          {
            for (vector<PeptideIdentification>::const_iterator pit =
                   citer->getPeptideIdentifications().begin(); pit !=
                 citer->getPeptideIdentifications().end(); ++pit)
            {
              writePeptideId(output, *pit, "PEPTIDE", false, false, false, peptide_id_meta_keys, peptide_hit_meta_keys);
            }
          }
        }
        outstr.close();
      }
      else if (in_type == FileTypes::CONSENSUSXML)
      {
        String consensus_centroids = getStringOption_("consensus:centroids");
        String consensus_elements = getStringOption_("consensus:elements");
        String consensus_features = getStringOption_("consensus:features");
        String sorting_method = getStringOption_("consensus:sorting_method");
        bool sort_by_maps = getFlag_("consensus:sort_by_maps");
        bool sort_by_size = getFlag_("consensus:sort_by_size");

        ConsensusMap consensus_map;
        ConsensusXMLFile consensus_xml_file;

        consensus_xml_file.load(in, consensus_map);

        // extract common id and hit meta values
        StringList peptide_id_meta_keys;
        StringList peptide_hit_meta_keys;

        vector<PeptideIdentification> pids;
        if (add_id_metavalues >= 0 || add_hit_metavalues >= 0)
        {
          const vector<PeptideIdentification>& uapids = consensus_map.getUnassignedPeptideIdentifications();
          pids.insert(pids.end(), uapids.begin(), uapids.end());
          for (ConsensusMap::const_iterator cmit = consensus_map.begin(); cmit != consensus_map.end(); ++cmit)
          {
              const vector<PeptideIdentification>& cpids = cmit->getPeptideIdentifications();
              pids.insert(pids.end(), cpids.begin(), cpids.end());
          }
          if (add_id_metavalues >= 0)
          {
            peptide_id_meta_keys = MetaInfoInterfaceUtils::findCommonMetaKeys<vector<PeptideIdentification>, StringList>(pids.begin(), pids.end(), add_id_metavalues);
              // currently there is some hardcoded logic to create extra columns for these meta values so remove them to prevent duplication
              peptide_id_meta_keys.erase(std::remove(peptide_id_meta_keys.begin(), peptide_id_meta_keys.end(), "predicted_RT"), peptide_id_meta_keys.end());
              peptide_id_meta_keys.erase(std::remove(peptide_id_meta_keys.begin(), peptide_id_meta_keys.end(), "predicted_RT_first_dim"), peptide_id_meta_keys.end());
              peptide_id_meta_keys.erase(std::remove(peptide_id_meta_keys.begin(), peptide_id_meta_keys.end(), "first_dim_rt"), peptide_id_meta_keys.end());
              peptide_id_meta_keys.erase(std::remove(peptide_id_meta_keys.begin(), peptide_id_meta_keys.end(), "predicted_PT"), peptide_id_meta_keys.end());
          }
          if (add_hit_metavalues >= 0)
          {
            vector<PeptideHit> temp_hits;
            for (Size i = 0; i != pids.size(); ++i)
            {
              const vector<PeptideHit>& hits = pids[i].getHits();
              temp_hits.insert(temp_hits.end(), hits.begin(), hits.end());
            }

            // siehe oben / analog machen
              peptide_hit_meta_keys = MetaInfoInterfaceUtils::findCommonMetaKeys<vector<PeptideHit>, StringList>(temp_hits.begin(), temp_hits.end(), add_hit_metavalues);
          }
        }

        if (sorting_method == "none")
        {
          // don't sort in this case
        }
        else if (sorting_method == "RT") consensus_map.sortByRT();
        else if (sorting_method == "MZ") consensus_map.sortByMZ();
        else if (sorting_method == "RT_then_MZ")
        {
          consensus_map.sortByPosition();
        }
        else if (sorting_method == "intensity")
        {
          consensus_map.sortByIntensity();
        }
        else if (sorting_method == "quality_decreasing")
        {
          consensus_map.sortByQuality(true);
        }
        else if (sorting_method == "quality_increasing")
        {
          consensus_map.sortByQuality(false);
        }

        if (sort_by_maps) consensus_map.sortByMaps();

        if (sort_by_size) consensus_map.sortBySize();

        String date_time_now = DateTime::now().get();

        // -------------------------------------------------------------------

        if (!consensus_centroids.empty())
        {
          std::ofstream consensus_centroids_file(consensus_centroids.c_str());
          if (!consensus_centroids_file)
          {
            throw Exception::UnableToCreateFile(__FILE__, __LINE__,
                                                OPENMS_PRETTY_FUNCTION,
                                                consensus_centroids);
          }

          SVOutStream output(consensus_centroids_file, sep, replacement,
                             quoting_method);

          writeConsensusHeader(output, "Centroids of consensus features", in,
                               date_time_now);
          writeFeatureHeader(output);
          output << nl;

          for (ConsensusMap::const_iterator cmit = consensus_map.begin();
               cmit != consensus_map.end(); ++cmit)
          {
            output << *cmit << nl;
          }
          consensus_centroids_file.close();
        }

        // -------------------------------------------------------------------

        if (!consensus_elements.empty())
        {
          std::ofstream consensus_elements_file(consensus_elements.c_str());
          if (!consensus_elements_file)
          {
            throw Exception::UnableToCreateFile(__FILE__, __LINE__,
                                                OPENMS_PRETTY_FUNCTION,
                                                consensus_elements);
          }

          SVOutStream output(consensus_elements_file, sep, replacement,
                             quoting_method);

          output.modifyStrings(false);
          writeConsensusHeader(output, "Elements of consensus features", in,
                               date_time_now);
          output << "#HL";
          writeFeatureHeader(output, "", false, false);
          writeFeatureHeader(output, "_cf", true, false);
          output << nl;
          output.modifyStrings(true);

          for (ConsensusMap::const_iterator cmit = consensus_map.begin();
               cmit != consensus_map.end(); ++cmit)
          {
            for (ConsensusFeature::const_iterator cfit = cmit->begin();
                 cfit != cmit->end(); ++cfit)
            {
              output << "H" << *cfit << *cmit << nl;
            }
            // We repeat the first feature handle at the end of the list.
            // This way you can generate closed line drawings
            // See Gnuplot set datafile commentschars
            output << "L" << *cmit->begin() << *cmit << nl;
          }
          consensus_elements_file.close();
        }

        // -------------------------------------------------------------------

        if (!consensus_features.empty())
        {
          std::ofstream consensus_features_file(consensus_features.c_str());
          if (!consensus_features_file)
          {
            throw Exception::UnableToCreateFile(__FILE__, __LINE__,
                                                OPENMS_PRETTY_FUNCTION,
                                                consensus_features);
          }

          SVOutStream output(consensus_features_file, sep, replacement,
                             quoting_method);

          std::map<Size, Size> map_id_to_map_num;
          std::vector<Size> map_num_to_map_id;
          FeatureHandle feature_handle_NaN;
          feature_handle_NaN.setRT(
            std::numeric_limits<FeatureHandle::CoordinateType>::quiet_NaN());
          feature_handle_NaN.setMZ(
            std::numeric_limits<FeatureHandle::CoordinateType>::quiet_NaN());
          feature_handle_NaN.setIntensity(
            std::numeric_limits<FeatureHandle::IntensityType>::quiet_NaN());
          // feature_handle_NaN.setCharge(std::numeric_limits<Int>::max());

          for (ConsensusMap::FileDescriptions::const_iterator fdit =
                 consensus_map.getFileDescriptions().begin();
               fdit != consensus_map.getFileDescriptions().end(); ++fdit)
          {
            map_id_to_map_num[fdit->first] = map_num_to_map_id.size();
            map_num_to_map_id.push_back(fdit->first);
          }

          map<String, Size> prot_runs;
          Size max_prot_run = 0;
          StringList comments;
          if (!no_ids)
          {
            String pep_line = "Protein identification runs associated with peptide/protein columns below: ";
            for (vector<ProteinIdentification>::const_iterator prot_it =
                   consensus_map.getProteinIdentifications().begin();
                 prot_it != consensus_map.getProteinIdentifications().end();
                 ++prot_it, ++max_prot_run)
            {
              String run_id = prot_it->getIdentifier();
              // add to comment:
              if (max_prot_run > 0)
              {
                pep_line += ", ";
              }
              pep_line += String(max_prot_run) + ": '" + run_id + "'";

              map<String, Size>::iterator pos = prot_runs.find(run_id);
              if (pos != prot_runs.end())
              {
                cerr << "Warning while exporting '" << in
                     << "': protein identification run ID '" << run_id
                     << "' occurs more than once" << endl;
              }
              else prot_runs[run_id] = max_prot_run;
            }
            if (max_prot_run > 0)
            {
              --max_prot_run; // increased beyond max. at end of for-loop
            }
            comments.push_back(pep_line);
          }

          writeConsensusHeader(output, "Consensus features", in,
                               date_time_now, comments);
          writeFeatureHeader(output, "_cf");
          output.modifyStrings(false);
          for (Size fhindex = 0; fhindex < map_num_to_map_id.size();
               ++fhindex)
          {
            Size map_id = map_num_to_map_id[fhindex];
            writeFeatureHeader(output, "_" + String(map_id), false, false);
          }
          if (!no_ids)
          {
            for (Size i = 0; i <= max_prot_run; ++i)
            {
              output << "peptide_" + String(i)
                     << "n_diff_peptides_" + String(i)
                     << "protein_" + String(i)
                     << "n_diff_proteins_" + String(i);
            }
          }
          output << nl;
          output.modifyStrings(true);

          for (ConsensusMap::const_iterator cmit = consensus_map.begin();
               cmit != consensus_map.end(); ++cmit)
          {
            output << *cmit;
            std::vector<FeatureHandle> feature_handles(map_num_to_map_id.size(),
                                                       feature_handle_NaN);
            for (ConsensusFeature::const_iterator cfit = cmit->begin();
                 cfit != cmit->end(); ++cfit)
            {
              feature_handles[map_id_to_map_num[cfit->getMapIndex()]] = *cfit;
            }
            for (Size fhindex = 0; fhindex < feature_handles.size();
                 ++fhindex)
            {
              output << feature_handles[fhindex];
            }
            if (!no_ids)
            {
              vector<set<String> > peptides_by_source(max_prot_run + 1),
              proteins_by_source(max_prot_run + 1);
              for (vector<PeptideIdentification>::const_iterator pep_it =
                     cmit->getPeptideIdentifications().begin(); pep_it !=
                   cmit->getPeptideIdentifications().end(); ++pep_it)
              {
                Size index = prot_runs[pep_it->getIdentifier()];
                for (vector<PeptideHit>::const_iterator hit_it = pep_it->
                                                                 getHits().begin(); hit_it != pep_it->getHits().end();
                     ++hit_it)
                {
                  peptides_by_source[index].insert(hit_it->getSequence().toString());
                  set<String> protein_accessions = hit_it->extractProteinAccessionsSet();
                  proteins_by_source[index].insert(protein_accessions.begin(), protein_accessions.end());
                }
              }
              vector<set<String> >::iterator pep_it = peptides_by_source.begin(), prot_it = proteins_by_source.begin();
              for (; pep_it != peptides_by_source.end(); ++pep_it, ++prot_it)
              {
                StringList seqs(vector<String>(pep_it->begin(),
                                               pep_it->end())),
                accs(vector<String>(prot_it->begin(), prot_it->end()));
                for (StringList::iterator acc_it = accs.begin();
                     acc_it != accs.end(); ++acc_it)
                {
                  acc_it->substitute('/', '_');
                }
                output << ListUtils::concatenate(seqs, "/") << seqs.size()
                       << ListUtils::concatenate(accs, "/") << accs.size();
              }
            }
            output << nl;
          }
          consensus_features_file.close();
        }

        // -------------------------------------------------------------------

        if (!out.empty())
        {
          std::ofstream outstr(out.c_str());
          if (!outstr)
          {
            throw Exception::UnableToCreateFile(__FILE__, __LINE__,
                                                OPENMS_PRETTY_FUNCTION, out);
          }

          SVOutStream output(outstr, sep, replacement, quoting_method);
          output.modifyStrings(false);
          writeConsensusHeader(output, "Consensus features", in,
                               date_time_now);

          std::map<Size, Size> map_id_to_map_num;
          std::vector<Size> map_num_to_map_id;
          FeatureHandle feature_handle_NaN;
          feature_handle_NaN.setRT(std::numeric_limits<
                                     FeatureHandle::CoordinateType>::quiet_NaN());
          feature_handle_NaN.setMZ(std::numeric_limits<
                                     FeatureHandle::CoordinateType>::quiet_NaN());
          feature_handle_NaN.setIntensity(std::numeric_limits<FeatureHandle::IntensityType>::quiet_NaN());
          feature_handle_NaN.setWidth(std::numeric_limits<
                                        FeatureHandle::WidthType>::quiet_NaN());
          feature_handle_NaN.setCharge(0); // just to be sure...
          // alternative?:
          // feature_handle_NaN.setCharge(std::numeric_limits<Int>::max());

          // It's hard to predict which meta keys will be used in file
          // descriptions. So we assemble a list each time. Represent keys
          // by String, not UInt, for implicit sorting.
          std::set<String> all_file_desc_meta_keys;
          std::vector<UInt> tmp_meta_keys;
          for (ConsensusMap::FileDescriptions::const_iterator fdit =
                 consensus_map.getFileDescriptions().begin();
               fdit != consensus_map.getFileDescriptions().end(); ++fdit)
          {
            map_id_to_map_num[fdit->first] = map_num_to_map_id.size();
            map_num_to_map_id.push_back(fdit->first);
            fdit->second.getKeys(tmp_meta_keys);
            for (std::vector<UInt>::const_iterator kit =
                   tmp_meta_keys.begin(); kit != tmp_meta_keys.end(); ++kit)
            {
              all_file_desc_meta_keys.insert(
                MetaInfoInterface::metaRegistry().getName(*kit));
            }
          }

          // headers (same order as the content of the output):
          output << "#MAP" << "id" << "filename" << "label" << "size";
          for (std::set<String>::const_iterator kit =
                 all_file_desc_meta_keys.begin(); kit !=
               all_file_desc_meta_keys.end(); ++kit)
          {
            output << *kit;
          }
          output << nl;
          if (!no_ids)
          {
            writeRunHeader(output);
            writeProteinHeader(output);
            writePeptideHeader(output, "UNASSIGNEDPEPTIDE");
            writeMetaValuesHeader(output, peptide_id_meta_keys);
            writeMetaValuesHeader(output, peptide_hit_meta_keys);
            output << nl;
          }
          output << "#CONSENSUS";
          writeFeatureHeader(output, "_cf", true, false);
          for (Size fhindex = 0; fhindex < map_num_to_map_id.size();
               ++fhindex)
          {
            Size map_id = map_num_to_map_id[fhindex];
            writeFeatureHeader(output, "_" + String(map_id), false, false);
          }
          output << nl;
          if (!no_ids)
          {
            writePeptideHeader(output, "PEPTIDE");
      writeMetaValuesHeader(output, peptide_id_meta_keys);
      writeMetaValuesHeader(output, peptide_hit_meta_keys);
            output << nl;
          }
          output.modifyStrings(true);

          // list of maps (intentionally at the beginning, contrary to order in consensusXML)
          for (ConsensusMap::FileDescriptions::const_iterator fdit =
                 consensus_map.getFileDescriptions().begin(); fdit !=
               consensus_map.getFileDescriptions().end(); ++fdit)
          {
            output << "MAP" << fdit->first << fdit->second.filename
                   << fdit->second.label << fdit->second.size;
            for (std::set<String>::const_iterator kit =
                   all_file_desc_meta_keys.begin(); kit !=
                 all_file_desc_meta_keys.end(); ++kit)
            {
              if (fdit->second.metaValueExists(*kit))
              {
                output << fdit->second.getMetaValue(*kit);
              }
              else output << "";
            }
            output << nl;
          }

          // proteins and unassigned peptides
          if (!no_ids) // proteins
          {
            for (vector<ProteinIdentification>::const_iterator it =
                   consensus_map.getProteinIdentifications().begin(); it !=
                 consensus_map.getProteinIdentifications().end(); ++it)
            {
              writeProteinId(output, *it);
            }

            // unassigned peptides
            for (vector<PeptideIdentification>::const_iterator pit = consensus_map.getUnassignedPeptideIdentifications().begin(); pit != consensus_map.getUnassignedPeptideIdentifications().end(); ++pit)
            {
              writePeptideId(output, *pit, "UNASSIGNEDPEPTIDE", false, false, false, peptide_id_meta_keys, peptide_hit_meta_keys);
              // first_dim_... stuff not supported for now
            }
          }

          // consensus features (incl. peptide annotations):
          for (ConsensusMap::const_iterator cmit = consensus_map.begin();
               cmit != consensus_map.end(); ++cmit)
          {
            std::vector<FeatureHandle> feature_handles(map_num_to_map_id.size(),
                                                       feature_handle_NaN);
            output << "CONSENSUS" << *cmit;
            for (ConsensusFeature::const_iterator cfit = cmit->begin();
                 cfit != cmit->end(); ++cfit)
            {
              feature_handles[map_id_to_map_num[cfit->getMapIndex()]] = *cfit;
            }
            for (Size fhindex = 0; fhindex < feature_handles.size(); ++fhindex)
            {
              output << feature_handles[fhindex];
            }
            output << nl;

            // peptide ids
            if (!no_ids)
            {
              for (vector<PeptideIdentification>::const_iterator pit =
                     cmit->getPeptideIdentifications().begin(); pit !=
                   cmit->getPeptideIdentifications().end(); ++pit)
              {
                writePeptideId(output, *pit, "PEPTIDE", false, false, false, peptide_id_meta_keys, peptide_hit_meta_keys);
              }
            }
          }
        }
        return EXECUTION_OK;
      }
      else if (in_type == FileTypes::IDXML)
      {
        vector<ProteinIdentification> prot_ids;
        vector<PeptideIdentification> pep_ids;
        String document_id;
        IdXMLFile().load(in, prot_ids, pep_ids, document_id);
        StringList peptide_id_meta_keys;
        StringList peptide_hit_meta_keys;

        if (add_id_metavalues >= 0) 
        {
          peptide_id_meta_keys = MetaInfoInterfaceUtils::findCommonMetaKeys<vector<PeptideIdentification>, StringList>(pep_ids.begin(), pep_ids.end(), add_id_metavalues);
          // currently there is some hardcoded logic to create extra columns for these meta values so remove them to prevent duplication 
          peptide_id_meta_keys.erase(std::remove(peptide_id_meta_keys.begin(), peptide_id_meta_keys.end(), "predicted_RT"), peptide_id_meta_keys.end());
          peptide_id_meta_keys.erase(std::remove(peptide_id_meta_keys.begin(), peptide_id_meta_keys.end(), "predicted_RT_first_dim"), peptide_id_meta_keys.end());
          peptide_id_meta_keys.erase(std::remove(peptide_id_meta_keys.begin(), peptide_id_meta_keys.end(), "first_dim_rt"), peptide_id_meta_keys.end());
          peptide_id_meta_keys.erase(std::remove(peptide_id_meta_keys.begin(), peptide_id_meta_keys.end(), "predicted_PT"), peptide_id_meta_keys.end());
        }

        if (add_hit_metavalues >= 0)
        {
          vector<PeptideHit> temp_hits;
          for (Size i = 0; i != pep_ids.size(); ++i)
          {
            const vector<PeptideHit>& hits = pep_ids[i].getHits();
            temp_hits.insert(temp_hits.end(), hits.begin(), hits.end());  
          }

          peptide_hit_meta_keys = MetaInfoInterfaceUtils::findCommonMetaKeys<vector<PeptideHit>, StringList>(temp_hits.begin(), temp_hits.end(), add_hit_metavalues);
        }

        ofstream txt_out(out.c_str());
        SVOutStream output(txt_out, sep, replacement, quoting_method);

        bool proteins_only = getFlag_("id:proteins_only");
        bool peptides_only = getFlag_("id:peptides_only");
        if (proteins_only && peptides_only)
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "'id:proteins_only' and 'id:peptides_only' cannot be used together");
        }

        String what = peptides_only ? "" : "PEPTIDE";
        if (!peptides_only)
        {
          writeRunHeader(output);
          writeProteinHeader(output);
        }
        if (!proteins_only)
        {
          writePeptideHeader(output, what, true, true, first_dim_rt);
          writeMetaValuesHeader(output, peptide_id_meta_keys);
          writeMetaValuesHeader(output, peptide_hit_meta_keys);
          output << nl;
        }

        for (vector<ProteinIdentification>::const_iterator it =
               prot_ids.begin(); it != prot_ids.end(); ++it)
        {
          String actual_id = it->getIdentifier();

          if (!peptides_only) writeProteinId(output, *it);

          if (!proteins_only)
          {
            // slight improvement on big idXML files with many different runs:
            // index the identifiers and peptide ids to avoid running over
            // them again and again (TODO)
            for (vector<PeptideIdentification>::const_iterator pit =
                   pep_ids.begin(); pit != pep_ids.end(); ++pit)
            {
              if (pit->getIdentifier() == actual_id)
              {
                writePeptideId(output, *pit, what, true, true, first_dim_rt, peptide_id_meta_keys, peptide_hit_meta_keys);
              }
            }
          }
        }

        txt_out.close();
      }
      else if (in_type == FileTypes::MZML)
      {
        PeakMap exp;
        FileHandler().loadExperiment(in, exp, FileTypes::MZML, ProgressLogger::NONE, false, false);

        if (exp.getSpectra().empty() && exp.getChromatograms().empty())
        {
          writeLog_("File does not contain spectra or chromatograms.");
          return INCOMPATIBLE_INPUT_DATA;
        }

        ofstream outstr(out.c_str());
        SVOutStream output(outstr, sep, replacement, quoting_method);
        output.modifyStrings(false);

        {
          if (exp.getSpectra().empty())
          {
            writeLog_("File does not contain spectra. No output for spectra generated!");
          }

          Size output_count(0);

          output << "#MS" << "level" << "rt" << "mz" << "charge" << "peaks" << "index" << "name" << nl;
          for (PeakMap::const_iterator it = exp.getSpectra().begin(); it != exp.getSpectra().end(); ++it)
          {
            int index = (it - exp.getSpectra().begin());
            String name = it->getName();
            if (it->getMSLevel() == 1)
            {
              ++output_count;
              output << "MS" << it->getMSLevel() << it->getRT() << "" << "" << it->size() << index << name << nl;
            }
            else if (it->getMSLevel() == 2)
            {
              double precursor_mz = -1;
              int precursor_charge = -1;

              if (!it->getPrecursors().empty())
              {
                precursor_mz = it->getPrecursors()[0].getMZ();
                precursor_charge = it->getPrecursors()[0].getCharge();
              }

              ++output_count;
              output << "MS" << it->getMSLevel() << it->getRT() << precursor_mz << precursor_charge << it->size() << index << name << nl;
            }
          }

          if (output_count != 0)
          {
            writeLog_("Exported " + String(output_count) + " spectra!");
          }
        }

        {
          if (exp.getChromatograms().empty())
          {
            writeLog_("File does not contain chromatograms. No output for chromatograms generated!");
          }

          Size output_count(0);
          Size unsupported_chromatogram_count(0);

          for (vector<MSChromatogram >::const_iterator it = exp.getChromatograms().begin(); it != exp.getChromatograms().end(); ++it)
          {
            if (it->getChromatogramType() == ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
            {
              ++output_count;
              output << "MRM Q1=" << it->getPrecursor().getMZ() << " Q3=" << it->getProduct().getMZ() << nl;
              for (MSChromatogram::ConstIterator cit = it->begin(); cit != it->end(); ++cit)
              {
                output << cit->getRT() << " " << cit->getIntensity() << nl;
              }
              output << nl;
            }
            else
            {
              ++unsupported_chromatogram_count;
            }
          }

          if (output_count != 0)
          {
            writeLog_("Exported " + String(output_count) + " SRM spectra!");
          }

          if (unsupported_chromatogram_count != 0)
          {
            writeLog_("Ignored " + String(unsupported_chromatogram_count) + " chromatograms not supported by TextExporter!");
          }
        }

        output << nl;
        outstr.close();
      }

      return EXECUTION_OK;
    }

  };
}


int main(int argc, const char** argv)
{
  TOPPTextExporter t;
  return t.main(argc, argv);
}

/// @endcond
