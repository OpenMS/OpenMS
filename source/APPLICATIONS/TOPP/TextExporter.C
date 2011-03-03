// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow, Marc Sturm, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/SVOutStream.h>

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

	 With the @p minimal flag, only the @p rt, @p mz, and @p intensity columns of @p FEATURE lines are written.

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
	 - one line per consensus feature (suitable for processing with e.g. @link http://www.r-project.org R@endlink)
	 - columns: same as for a @p CONSENSUS line above, followed by additional columns for identification data
	 - additional columns: @p peptide_N0, @p n_diff_peptides_N0, @p protein_N0, @p n_diff_proteins_N0, @p peptide_N1, ...
	 - @p "..._Ni" columns refer to the identification run with index "Ni", @p n_diff_... stands for "number of different ..."; different peptides/proteins in one column are separated by "/"

	 With the @p no_ids flag, the additional columns are not included.

	 <B>idXML input:</B>
	 - first column: @p RUN / @p PROTEIN / @p PEPTIDE (indicator for the type of data in the current row)
	 - see above for the formats of @p RUN, @p PROTEIN, @p PEPTIDE lines
	 - additional column for @p PEPTIDE lines: @p predicted_rt

	 With the @p proteins_only flag, only @p RUN and @p PROTEIN lines are written.

	 With the @p peptides_only flag, only @p PEPTIDE lines (without the @p PEPTIDE indicator) are written.

	 With the @p first_dim_rt flag, the additional columns @p rt_first_dim and @p predicted_rt_first_dim are included for @p PEPTIDE lines.

	 
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_TextExporter.cli
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
		StringList elements = StringList::create("#rt,mz,intensity,charge,width");
		if (!comment) elements[0] = "rt";
		if (incl_quality) elements << "quality";
		bool old = out.modifyStrings(false);
		for (StringList::iterator it = elements.begin(); it != elements.end(); ++it)
		{
			out << *it + suffix;
		}
		out.modifyStrings(old);
	}

	// write the header for exporting consensusXML
	void writeConsensusHeader(SVOutStream& out, const String& what, 
														const String& infile,	const String& now,
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
				<< "date_time" << "search_engine_version" << "parameters" << endl;
		out.modifyStrings(old);
	}

	// write the header for protein data
	void writeProteinHeader(SVOutStream& out)
	{
		bool old = out.modifyStrings(false);
		out << "#PROTEIN" << "score" << "rank" << "accession" << "coverage" << "sequence" << endl;
		out.modifyStrings(old);
	}

	// stream output operator for a ProteinHit
	SVOutStream& operator<<(SVOutStream& out, const ProteinHit& hit)
	{
    out << hit.getScore() << hit.getRank() << hit.getAccession()
        << hit.getCoverage() << hit.getSequence();
		return out;
	}

	// stream output operator for SearchParameters
	SVOutStream& operator<<(SVOutStream& out,
													const ProteinIdentification::SearchParameters sp)
	{
		String param_line = "db=" + sp.db + ", db_version=" +	sp.db_version +
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
    param_line += ProteinIdentification::NamesOfDigestionEnzyme[sp.enzyme];
		param_line += ", missed_cleavages=" + String(sp.missed_cleavages) +
			", peak_mass_tolerance=" + String(sp.peak_mass_tolerance) +
			", precursor_mass_tolerance=" + String(sp.precursor_tolerance);
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
		out << sp << endl;
		for (vector<ProteinHit>::const_iterator hit_it = pid.getHits().begin();
				 hit_it != pid.getHits().end(); ++hit_it)
		{
			out << "PROTEIN" << *hit_it << endl;
		}
	}

	// write the header for peptide data
	void writePeptideHeader(SVOutStream& out, const String& what = "PEPTIDE",
													bool incl_pred_rt = false,
													bool incl_first_dim = false)
	{
		bool old = out.modifyStrings(false);
		if (what.empty()) out << "#rt";
		else out << "#" + what << "rt";
		out << "mz" << "score" << "rank" << "sequence" << "charge" << "aa_before" 
				<< "aa_after" << "score_type" << "search_identifier" << "accessions";
		if (incl_pred_rt) out << "predicted_rt";
		if (incl_first_dim) out << "rt_first_dim" << "predicted_rt_first_dim";
		out << endl;
		out.modifyStrings(old);
	}

	// stream output operator for a PeptideHit
	SVOutStream& operator<<(SVOutStream& out, const PeptideHit& hit)
	{
		out << hit.getScore() << hit.getRank() << hit.getSequence()
				<< hit.getCharge() << hit.getAABefore() << hit.getAAAfter();
		return out;
	}

  // write a protein identification to the output stream
	void writePeptideId(SVOutStream& out, const PeptideIdentification& pid,
											const String& what = "PEPTIDE", bool incl_pred_rt = false,
											bool incl_first_dim = false)
	{
		for (vector<PeptideHit>::const_iterator hit_it = pid.getHits().begin();
				 hit_it != pid.getHits().end(); ++hit_it)
		{
			if (!what.empty()) out << what;
			if (pid.metaValueExists("RT")) out << (DoubleReal) pid.getMetaValue("RT");
			else out << "-1";
			if (pid.metaValueExists("MZ")) out << (DoubleReal) pid.getMetaValue("MZ");
			else out << "-1";
			out << *hit_it << pid.getScoreType() << pid.getIdentifier();
			String accessions;
			for (vector<String>::const_iterator acc_it =
						 hit_it->getProteinAccessions().begin(); acc_it !=
						 hit_it->getProteinAccessions().end(); ++acc_it)
			{
				if (acc_it != hit_it->getProteinAccessions().begin())
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
			out << endl;
		}
	}


	class TOPPTextExporter : public TOPPBase
  {
	public:
		TOPPTextExporter() :
			TOPPBase("TextExporter", "Exports various XML formats to a text file.")
      {
      }

	protected:


		void registerOptionsAndFlags_()
      {
        registerInputFile_("in", "<file>", "", "Input file ");
        setValidFormats_("in", StringList::create(
          "featureXML,consensusXML,idXML,mzML"));
        registerOutputFile_("out", "<file>", "",
          "Output file (mandatory for featureXML and idXML)", false);
        registerStringOption_("separator", "<sep>", "\t", "The used separator character(s); if not set the 'tab' character is used", false);
				registerStringOption_("replacement", "<string>", "_", "Used to replace occurrences of the separator in strings before writing, if 'quoting' is 'none'", false);
				registerStringOption_("quoting", "<method>", "none", "Method for quoting of strings: 'none' for no quoting, 'double' for quoting with doubling of embedded quotes,\n'escape' for quoting with backslash-escaping of embedded quotes", false);
				setValidStrings_("quoting", StringList::create("none,double,escape"));
        registerFlag_("no_ids", "Suppresses output of identification data.");
        addEmptyLine_();

				addText_("Options for featureXML files:");
				registerFlag_("minimal", "Set this flag to write only three attributes: RT, m/z, and intensity.");
				addEmptyLine_();

        addText_("Options for idXML files:");
        registerFlag_("proteins_only",
          "Set this flag if you want only protein information from an idXML file");
        registerFlag_("peptides_only",
          "Set this flag if you want only peptide information from an idXML file");
        registerFlag_(
          "first_dim_rt",
          "If this flag is set the first_dim RT of the peptide hits will also be printed (if present).");
        addEmptyLine_();

        addText_("Options for consensusXML files:");
        registerOutputFile_("consensus_centroids", "<file>", "",
          "Output file for centroids of consensus features", false);
        registerOutputFile_("consensus_elements", "<file>", "",
          "Output file for elements of consensus features", false);
        registerOutputFile_("consensus_features", "<file>", "", "Output file for consensus features and contained elements from all maps (writes 'nan's if elements are missing)", false);
        registerStringOption_("sorting_method", "<method>", "none",
          "Sorting method", false);
        setValidStrings_("sorting_method", StringList::create("none,RT,MZ,RT_then_MZ,intensity,quality_decreasing,quality_increasing"));
        registerFlag_("sort_by_maps",
          "Apply a stable sort by the covered maps, lexicographically", false);
        registerFlag_("sort_by_size", "Apply a stable sort by decreasing size (i.e., the number of elements)", false);
        addText_("Sorting options can be combined.  The precedence is: sort_by_size, sort_by_maps, sorting_method");
      }


		ExitCodes main_( int, const char** )
      {
        //-------------------------------------------------------------
        // parameter handling
        //-------------------------------------------------------------
        String in = getStringOption_("in");
        String out = getStringOption_("out");
        bool no_ids = getFlag_("no_ids");
        bool first_dim_rt = getFlag_("first_dim_rt");

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
        writeDebug_(String("Input file type: ") + FileHandler::typeToName(
          in_type), 2);

        if ( in_type == FileTypes::UNKNOWN )
        {
          writeLog_("Error: Could not determine input file type!");
          return PARSE_ERROR;
        }

        if ( in_type == FileTypes::FEATUREXML )
        {
          //-------------------------------------------------------------
          // loading input
          //-------------------------------------------------------------

          FeatureMap<> feature_map;
          FeatureXMLFile f;
          f.load(in, feature_map);

          // compute protein coverage
          vector<ProteinIdentification> prot_ids = feature_map.getProteinIdentifications();
          vector<PeptideIdentification> pep_ids;
					// collect all peptide ids:
          for (Size i = 0; i < feature_map.size(); ++i)
          {
            vector<PeptideIdentification> pep_ids_bf = feature_map[i].getPeptideIdentifications();
            pep_ids.insert(pep_ids.end(), pep_ids_bf.begin(), pep_ids_bf.end());
          }
          pep_ids.insert(pep_ids.end(), feature_map.getUnassignedPeptideIdentifications().begin(), feature_map.getUnassignedPeptideIdentifications().end());

          try // might throw Exception::MissingInformation()
          {
            for (Size i = 0; i < prot_ids.size(); ++i) 
						{
							prot_ids[i].computeCoverage(pep_ids);
						}
          }
          catch (...) {}
          feature_map.setProteinIdentifications(prot_ids);

          // text output
          ofstream outstr(out.c_str());
					SVOutStream output(outstr, sep, replacement, quoting_method);

					bool minimal = getFlag_("minimal");
					no_ids |= minimal; // "minimal" implies "no_ids"

					// write header:
					output.modifyStrings(false);
					bool comment = true;
					if (!no_ids) 
					{
						writeRunHeader(output);
						writeProteinHeader(output);
						writePeptideHeader(output, "UNASSIGNEDPEPTIDE");
						output << "#FEATURE";
						comment = false;
					}
					if (minimal) output << "#rt" << "mz" << "intensity";
					else
					{
						writeFeatureHeader(output, "", true, comment);
						output << "rt_quality" << "mz_quality" << "rt_start" << "rt_end";
					}
					output << endl;
					if (!no_ids)
					{
						writePeptideHeader(output);
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
							writePeptideId(output, *pit, "UNASSIGNEDPEPTIDE");
						}
					}

          for (FeatureMap<>::const_iterator citer = feature_map.begin(); 
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
								output << citer->getConvexHulls().begin()->
									getBoundingBox().minX() << citer->getConvexHulls().begin()->
									getBoundingBox().maxX();
							}
							else 
							{
								output << "-1" << "-1";
							}
						}
            output << endl;

            //peptide ids
            if (!no_ids)
            {
              for (vector<PeptideIdentification>::const_iterator pit =
										 citer->getPeptideIdentifications().begin(); pit !=
										 citer->getPeptideIdentifications().end(); ++pit)
              {
								writePeptideId(output, *pit);
							}
						}
					}
					outstr.close();
				}

        else if (in_type == FileTypes::CONSENSUSXML)
        {
          String consensus_centroids = getStringOption_("consensus_centroids");
          String consensus_elements = getStringOption_("consensus_elements");
          String consensus_features = getStringOption_("consensus_features");
          String sorting_method = getStringOption_("sorting_method");
          bool sort_by_maps = getFlag_("sort_by_maps");
          bool sort_by_size = getFlag_("sort_by_size");

          ConsensusMap consensus_map;
          ConsensusXMLFile consensus_xml_file;

          consensus_xml_file.load(in, consensus_map);

          // compute protein coverage
          vector<ProteinIdentification> prot_ids = consensus_map.getProteinIdentifications();
          vector<PeptideIdentification> pep_ids;
          for (Size i=0;i<consensus_map.size();++i) // collect all peptide ids
          {
            vector<PeptideIdentification> pep_ids_bf = consensus_map[i].getPeptideIdentifications();
            pep_ids.insert(pep_ids.end(), pep_ids_bf.begin(), pep_ids_bf.end());
          }
          pep_ids.insert(pep_ids.end(), consensus_map.getUnassignedPeptideIdentifications().begin(), consensus_map.getUnassignedPeptideIdentifications().end());
          try
          { // might throw Exception::MissingInformation()
            for (Size i=0;i<prot_ids.size();++i) prot_ids[i].computeCoverage(pep_ids);
          }
          catch (...){}
          consensus_map.setProteinIdentifications(prot_ids);


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
                __PRETTY_FUNCTION__, consensus_centroids);
            }

						SVOutStream output(consensus_centroids_file, sep, replacement,
															 quoting_method);

						writeConsensusHeader(output, "Centroids of consensus features", in,
																 date_time_now);
						writeFeatureHeader(output);
						output << endl;

            for (ConsensusMap::const_iterator cmit = consensus_map.begin();
								 cmit != consensus_map.end(); ++cmit)
            {
              output << *cmit << endl;
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
                __PRETTY_FUNCTION__, consensus_elements);
            }

						SVOutStream output(consensus_elements_file, sep, replacement,
															 quoting_method);

						output.modifyStrings(false);
						writeConsensusHeader(output, "Elements of consensus features", in,
																 date_time_now);
						output << "#HL";
						writeFeatureHeader(output, "", false, false);
						writeFeatureHeader(output, "_cf", true, false);
						output << endl;
						output.modifyStrings(true);

            for (ConsensusMap::const_iterator cmit = consensus_map.begin();
								 cmit != consensus_map.end(); ++cmit)
            {
              for (ConsensusFeature::const_iterator cfit = cmit->begin();
									 cfit != cmit->end(); ++cfit)
              {
                output << "H" << *cfit << *cmit << endl;
              }
              // We repeat the first feature handle at the end of the list.
              // This way you can generate closed line drawings
              // See Gnuplot set datafile commentschars
              output << "L" << *cmit->begin() << *cmit << endl;
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
                __PRETTY_FUNCTION__, consensus_features);
            }

						SVOutStream output(consensus_features_file, sep, replacement,
															 quoting_method);

						std::map<Size,Size> map_id_to_map_num;
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
								 fdit != consensus_map.getFileDescriptions().end(); ++fdit )
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
							if (max_prot_run>0)
							{
								--max_prot_run; // increased beyond max. at end of for-loop
							}
							comments << pep_line;
						}

						writeConsensusHeader(output, "Consensus features", in, 
																 date_time_now,	comments);
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
            output << endl;
						output.modifyStrings(true);

            for (ConsensusMap::const_iterator cmit = consensus_map.begin();
								 cmit != consensus_map.end(); ++cmit)
            {
              output << *cmit;
              std::vector<FeatureHandle> feature_handles(map_num_to_map_id.size(), feature_handle_NaN);
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
										peptides_by_source[index].insert(hit_it->getSequence().
																										 toString());
										proteins_by_source[index].insert(
											hit_it->getProteinAccessions().begin(),
											hit_it->getProteinAccessions().end());
									}
								}
								vector<set<String> >::iterator pep_it = peptides_by_source.
									begin(), prot_it = proteins_by_source.begin();
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
									output << seqs.concatenate("/") << seqs.size()
												 << accs.concatenate("/") << accs.size();
								}
							}
              output << endl;
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
                __PRETTY_FUNCTION__, out);
            }

						SVOutStream output(outstr, sep, replacement, quoting_method);
						output.modifyStrings(false);
						writeConsensusHeader(output, "Consensus features", in, 
																 date_time_now);

            std::map<Size,Size> map_id_to_map_num;
            std::vector<Size> map_num_to_map_id;
            FeatureHandle feature_handle_NaN;
            feature_handle_NaN.setRT(std::numeric_limits<
                FeatureHandle::CoordinateType>::quiet_NaN());
            feature_handle_NaN.setMZ(std::numeric_limits<
                FeatureHandle::CoordinateType>::quiet_NaN());
            feature_handle_NaN.setIntensity(std::numeric_limits<
                FeatureHandle::IntensityType>::quiet_NaN());
            feature_handle_NaN.setWidth(std::numeric_limits<
                FeatureHandle::WidthType>::quiet_NaN());
            feature_handle_NaN.setCharge(0); // just to be sure...
            // feature_handle_NaN.setCharge(std::numeric_limits<Int>::max()); // alternative ??

            // It's hard to predict which meta keys will be used in file
						// descriptions. So we assemble a list each time. Represent keys
						// by String, not UInt, for implicit sorting.
            std::set<String> all_file_desc_meta_keys;
            std::vector<UInt> tmp_meta_keys;
            for (ConsensusMap::FileDescriptions::const_iterator fdit =
								 consensus_map.getFileDescriptions().begin(); 
								 fdit != consensus_map.getFileDescriptions().end(); ++fdit )
            {
              map_id_to_map_num[fdit->first] = map_num_to_map_id.size();
              map_num_to_map_id.push_back(fdit->first);
              fdit->second.getKeys(tmp_meta_keys);
              for (std::vector<UInt>::const_iterator kit =
										 tmp_meta_keys.begin(); kit != tmp_meta_keys.end(); ++kit )
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
            output << endl;
            if (!no_ids)
            {
							writeRunHeader(output);
							writeProteinHeader(output);
							writePeptideHeader(output, "UNASSIGNEDPEPTIDE");
						}
						output << "#CONSENSUS";
						writeFeatureHeader(output, "_cf", true, false);
						for (Size fhindex = 0; fhindex < map_num_to_map_id.size();
								 ++fhindex)
						{
							Size map_id = map_num_to_map_id[fhindex];
							writeFeatureHeader(output, "_" + String(map_id), false, false);
						}
						output << endl;
						if (!no_ids) writePeptideHeader(output, "PEPTIDE");
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
              output << endl;
            }

            // proteins and unassigned peptides
            if (!no_ids)
            { // proteins
							for (vector<ProteinIdentification>::const_iterator it =
										 consensus_map.getProteinIdentifications().begin(); it !=
										 consensus_map.getProteinIdentifications().end(); ++it )
							{
								writeProteinId(output, *it);
							}

              // unassigned peptides
							for (vector<PeptideIdentification>::const_iterator pit = consensus_map.getUnassignedPeptideIdentifications().begin(); pit != consensus_map.getUnassignedPeptideIdentifications().end(); ++pit)
							{
								writePeptideId(output, *pit, "UNASSIGNEDPEPTIDE");
								// first_dim_... stuff not supported for now
							}
            }

						// consensus features (incl. peptide annotations):
            for (ConsensusMap::const_iterator cmit = consensus_map.begin();
								 cmit != consensus_map.end(); ++cmit)
            {
              std::vector<FeatureHandle> feature_handles(map_num_to_map_id.size(), feature_handle_NaN);
              output << "CONSENSUS" << *cmit;
              for (ConsensusFeature::const_iterator cfit = cmit->begin();
									 cfit != cmit->end(); ++cfit)
              {
                feature_handles[map_id_to_map_num[cfit->getMapIndex()]] = *cfit;
              }
              for (Size fhindex = 0; fhindex < feature_handles.size();
									 ++fhindex )
              {
								output << feature_handles[fhindex];
              }
              output << endl;

              // peptide ids
              if (!no_ids)
              {
                for (vector<PeptideIdentification>::const_iterator pit =
											 cmit->getPeptideIdentifications().begin(); pit !=
											 cmit->getPeptideIdentifications().end(); ++pit)
                {
									writePeptideId(output, *pit);
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
          try
          { // might throw Exception::MissingInformation()
            for (Size i = 0; i < prot_ids.size(); ++i)
						{
							prot_ids[i].computeCoverage(pep_ids);
						}
          }
          catch (...) {}

          ofstream txt_out(out.c_str());
					SVOutStream output(txt_out, sep, replacement, quoting_method);

					bool proteins_only = getFlag_("proteins_only");
					bool peptides_only = getFlag_("peptides_only");
					if (proteins_only && peptides_only)
					{
						throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "'proteins_only' and 'peptides_only' cannot be used together");
					}

					String what = peptides_only ? "" : "PEPTIDE";
					if (!peptides_only)
					{
						writeRunHeader(output);
						writeProteinHeader(output);
					}
					if (!proteins_only)
					{
						writePeptideHeader(output, what, true, first_dim_rt);
					}

          for (vector<ProteinIdentification>::const_iterator it =
								 prot_ids.begin(); it != prot_ids.end(); ++it)
          {
            String actual_id = it->getIdentifier();

            if (!peptides_only)	writeProteinId(output, *it);

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
									writePeptideId(output, *pit, what, true, first_dim_rt);
                }
              }
            }
          }

          txt_out.close();
        }

				else if (in_type == FileTypes::MZML)
				{
					PeakMap exp;
					FileHandler().loadExperiment(in, exp);

          if (exp.getChromatograms().size()==0)
          {
            writeLog_("File does not contain chromatograms. No output was generated!");
            return INCOMPATIBLE_INPUT_DATA;
          }

          Size output_count(0);

					ofstream outstr(out.c_str());
					SVOutStream output(outstr, sep, replacement, quoting_method);
					output.modifyStrings(false);
					for (vector<MSChromatogram<> >::const_iterator it = exp.getChromatograms().begin(); it != exp.getChromatograms().end(); ++it)
					{
						if (it->getChromatogramType() == ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
						{
              ++output_count;
							output << "MRM Q1=" << it->getPrecursor().getMZ() << " Q3=" << it->getProduct().getMZ() << endl;
							for (MSChromatogram<>::ConstIterator cit = it->begin(); cit != it->end(); ++cit)
							{
								output << cit->getRT() << " " << cit->getIntensity() << endl;
							}
							output << endl;
						}
					}
					outstr.close();

          writeLog_("Found " + String() + " SRM spectra!");
          if (output_count==0) writeLog_("No output was generated!!");
				}

        return EXECUTION_OK;
      }
  };
}


int
main( int argc, const char** argv )
{
  TOPPTextExporter t;
  return t.main(argc, argv);
}

/// @endcond
