// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Andreas Bertsch, Chris Bielow, Marc Sturm $
// $Authors: $
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

 @brief This application converts several %OpenMS XML formats
 (namely featureXML, consensusXML and idXML) to text files.

 The primary goal of this tool is to create a readable format
 for Excel and OpenOffice.

 <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_TextExporter.cli
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{

  // There is no standard how to print a nan (not-a-number) value.  So we do
  // this on our own.
  namespace
  {
    const char nan[] = "nan"; // that's what Linux GCC uses, and gnuplot understands.

    template < typename NumberT >
      std::ostream &
      printValueOrNan( std::ostream & os, NumberT thing )
      {
        if ( !boost::math::isnan(thing) )
        {
          return os << thing;
        }
        else
        {
          return os << nan;
        }
      }
  }

  /// Wrapper class to implement printing of FeatureHandle
  struct FeatureHandlePrinter
  {
      FeatureHandlePrinter( const FeatureHandle &rhs ) :
        ref_(rhs)
      {
      }
      const FeatureHandle &ref_;
  };

  /// Output operator for a FeatureHandlePrinter.
  std::ostream &
  operator <<( std::ostream& os, const FeatureHandlePrinter& rhs )
  {
    const Size exponent_extra_digits = 6;
    const Size charge_digits = 5;
    const unsigned prec_save = os.precision();
    os << std::setprecision(writtenDigits<> (FeatureHandle::CoordinateType()))
        << std::setw(writtenDigits(FeatureHandle::CoordinateType())
            + exponent_extra_digits);
    printValueOrNan(os, rhs.ref_.getRT());
    os << ' ' << std::setw(writtenDigits(FeatureHandle::CoordinateType())
        + exponent_extra_digits);
    printValueOrNan(os, rhs.ref_.getMZ());
    os << ' ' << std::setprecision(writtenDigits<> (
      FeatureHandle::IntensityType())) << std::setw(writtenDigits(
      FeatureHandle::IntensityType()) + exponent_extra_digits);
    printValueOrNan(os, rhs.ref_.getIntensity());
    os << ' ' << std::setw(charge_digits) << rhs.ref_.getCharge()
        << std::setprecision(prec_save);
    return os;
  }

  /// Wrapper class to implement printing of ConsensusFeature
  struct ConsensusFeaturePrinter
  {
      ConsensusFeaturePrinter( const ConsensusFeature &rhs ) :
        ref_(rhs)
      {
      }
      const ConsensusFeature &ref_;
  };

  /// Output operator for a ConsensusFeaturePrinter.
  std::ostream &
  operator <<( std::ostream& os, const ConsensusFeaturePrinter& rhs )
  {
    const Size exponent_extra_digits = 6;
    const Size charge_digits = 5;
    const unsigned prec_save = os.precision();
    os << std::setprecision(writtenDigits(FeatureHandle::CoordinateType()))
        << std::setw(writtenDigits(FeatureHandle::CoordinateType())
            + exponent_extra_digits);
    printValueOrNan(os, rhs.ref_.getRT());
    os << ' ' << std::setw(writtenDigits(FeatureHandle::CoordinateType())
        + exponent_extra_digits);
    printValueOrNan(os, rhs.ref_.getMZ());
    os << ' ' << std::setprecision(
      writtenDigits(FeatureHandle::IntensityType())) << std::setw(
      writtenDigits(FeatureHandle::IntensityType()) + exponent_extra_digits);
    printValueOrNan(os, rhs.ref_.getIntensity());
    os << ' ' << std::setw(charge_digits) << rhs.ref_.getCharge()
        << std::setprecision(prec_save);
    return os;
  }

  class TOPPTextExporter : public TOPPBase
  {
    public:
      TOPPTextExporter() :
        TOPPBase("TextExporter", "Exports various XML formats to a text file.")
      {
      }

    protected:

      void
      registerOptionsAndFlags_()
      {
        registerInputFile_("in", "<file>", "", "Input file ");
        setValidFormats_("in", StringList::create(
          "featureXML,consensusXML,idXML"));
        registerOutputFile_("out", "<file>", "",
          "Output file. Mandatory for FeatureXML and IdXML.", false);
        registerStringOption_(
          "separator",
          "<sep>",
          "",
          "The used separator characters. If unset the 'tab' character is used.",
          false);
        registerFlag_("no_ids",
          "Suppresses output of identification data.");
        addEmptyLine_();

        addText_("Options for IdXML files:");
        registerFlag_("proteins_only",
          "Set this flag if you want only protein information from an idXML file");
        registerFlag_("peptides_only",
          "Set this flag if you want only peptide information from an idXML file");
        registerFlag_(
          "first_dim_rt",
          "If this flag is set the first_dim RT of the peptide hits will also be printed (if present).");
        addEmptyLine_();

        addText_("Options for ConsensusXML files:");
        registerOutputFile_("consensus_centroids", "<file>", "",
          "Centroids of consensus features", false);
        registerOutputFile_("consensus_elements", "<file>", "",
          "Elements of consensus features", false);
        registerOutputFile_(
          "consensus_features",
          "<file>",
          "",
          "Consensus features and contained elements from all maps (writes 'nan's if element is missing)",
          false);
        addText_("Each of the consensus_... files is created as requested.");
        registerStringOption_("sorting_method", "<method>", "none",
          "Sorting method", false);
        setValidStrings_(
          "sorting_method",
          StringList::create(
            "none,RT,MZ,RT_then_MZ,intensity,quality_decreasing,quality_increasing"));
        registerFlag_("sort_by_maps",
          "Apply a stable sort by the covered maps, lexicographically", false);
        registerFlag_(
          "sort_by_size",
          "Apply a stable sort by decreasing size (i.e., the number of elements)",
          false);
        addText_(
          "Sorting options can be combined.  The precedence is: sort_by_size, sort_by_maps, sorting_method");
      }

      ExitCodes
      main_( int, const char** )
      {

        //-------------------------------------------------------------
        // parameter handling
        //-------------------------------------------------------------
        String in = getStringOption_("in");
        String out = getStringOption_("out");
        Size counter = 0;
        bool no_ids = getFlag_("no_ids");
        bool first_dim_rt = getFlag_("first_dim_rt");

        //separator
        String sep = getStringOption_("separator");
        if ( sep == "" ) sep = "\t";

        //input file type
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

          // text output
          ofstream outstr(out.c_str());

          // stores one feature per line
          if ( no_ids )
          {
            outstr << "#rt" << sep << "mz" << sep << "intensity" << sep
                << "charge" << sep << "overall_quality" << sep << "rt_quality"
                << sep << "mz_quality" << sep << "rt_start" << sep << "rt_end"
                << endl;
          }
          else
          {
            outstr << "#FEATURE" << sep << "rt" << sep << "mz" << sep
                << "intensity" << sep << "charge" << sep << "overall_quality"
                << sep << "rt_quality" << sep << "mz_quality" << sep
                << "rt_start" << sep << "rt_end" << endl;
            outstr << "#PEPTIDE" << sep << "rt" << sep << "mz" << sep
                << "score" << sep << "rank" << sep << "sequence" << sep
                << "charge" << sep << "AA_before" << sep << "AA_after" << sep
                << "score_type" << sep << "search_identifier" << endl;
          }
          for ( FeatureMap<>::const_iterator citer = feature_map.begin(); citer
              != feature_map.end(); ++citer )
          {
            if ( !no_ids ) outstr << "FEATURE" << sep;
            outstr << citer->getPosition()[0] << sep << citer->getPosition()[1]
                << sep << citer->getIntensity();
            outstr << sep << citer->getCharge();
            outstr << sep << citer->getOverallQuality();
            outstr << sep << citer->getQuality(0) << sep
                << citer->getQuality(1);

            if ( citer->getConvexHulls().size() > 0 )
            {
              outstr << sep
                  << citer->getConvexHulls().begin()->getBoundingBox().minX();
              outstr << sep
                  << citer->getConvexHulls().begin()->getBoundingBox().maxX();
            }
            else
            {
              outstr << sep << "-1";
              outstr << sep << "-1";
            }
            outstr << endl;

            //peptide ids
            if ( !no_ids )
            {
              for ( vector<PeptideIdentification>::const_iterator pit =
                  citer->getPeptideIdentifications().begin(); pit
                  != citer->getPeptideIdentifications().end(); ++pit )
              {
                for ( vector<PeptideHit>::const_iterator ppit =
                    pit->getHits().begin(); ppit != pit->getHits().end(); ++ppit )
                {
                  outstr << "PEPTIDE" << sep;
                  if ( pit->metaValueExists("RT") )
                  {
                    outstr << (DoubleReal) pit->getMetaValue("RT") << sep;
                  }
                  else
                  {
                    outstr << "-1" << sep;
                  }

                  if ( pit->metaValueExists("MZ") )
                  {
                    outstr << (DoubleReal) pit->getMetaValue("MZ") << sep;
                  }
                  else
                  {
                    outstr << "-1" << sep;
                  }
                  outstr << ppit->getScore() << sep << ppit->getRank() << sep
                      << ppit->getSequence() << sep << ppit->getCharge() << sep
                      << ppit->getAABefore() << sep << ppit->getAAAfter()
                      << sep << pit->getScoreType() << sep
                      << pit->getIdentifier() << endl;
                }
              }
            }
          }
          outstr.close();
        }
        else if ( in_type == FileTypes::CONSENSUSXML )
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

          if ( sorting_method == "none" )
          {
            // don't sort in this case
          }
          else if ( sorting_method == "RT" )
          {
            consensus_map.sortByRT();
          }
          else if ( sorting_method == "MZ" )
          {
            consensus_map.sortByMZ();
          }
          else if ( sorting_method == "RT_then_MZ" )
          {
            consensus_map.sortByPosition();
          }
          else if ( sorting_method == "intensity" )
          {
            consensus_map.sortByIntensity();
          }
          else if ( sorting_method == "quality_decreasing" )
          {
            consensus_map.sortByQuality(true);
          }
          else if ( sorting_method == "quality_increasing" )
          {
            consensus_map.sortByQuality(false);
          }

          if ( sort_by_maps )
          {
            consensus_map.sortByMaps();
          }

          if ( sort_by_size )
          {
            consensus_map.sortBySize();
          }

          String date_time_now = DateTime::now().get();

          // ----------------------------------------------------------------------

          if ( !consensus_centroids.empty() )
          {
            std::ofstream consensus_centroids_file(consensus_centroids.c_str());
            if ( !consensus_centroids_file )
            {
              throw Exception::UnableToCreateFile(__FILE__, __LINE__,
                __PRETTY_FUNCTION__, consensus_centroids);
            }

            consensus_centroids_file
                << "#  Centroids of consensus features extracted from " << in
                << " on " << date_time_now << std::endl
                << "# RT MZ Intensity Charge" << std::endl;
            for ( ConsensusMap::const_iterator cmit = consensus_map.begin(); cmit
                != consensus_map.end(); ++cmit )
            {
              consensus_centroids_file << ConsensusFeaturePrinter(*cmit)
                  << std::endl;
            }
            consensus_centroids_file.close();
          }

          // ----------------------------------------------------------------------

          if ( !consensus_elements.empty() )
          {
            std::ofstream consensus_elements_file(consensus_elements.c_str());
            if ( !consensus_elements_file )
            {
              throw Exception::UnableToCreateFile(__FILE__, __LINE__,
                __PRETTY_FUNCTION__, consensus_elements);
            }

            consensus_elements_file
                << "#  Elements of consensus features extracted from " << in
                << " on " << date_time_now << std::endl
                << "# RT MZ Intensity Charge" << std::endl;
            for ( ConsensusMap::const_iterator cmit = consensus_map.begin(); cmit
                != consensus_map.end(); ++cmit )
            {
              consensus_elements_file << std::endl;
              for ( ConsensusFeature::const_iterator cfit = cmit->begin(); cfit
                  != cmit->end(); ++cfit )
              {
                consensus_elements_file << "H " << FeatureHandlePrinter(*cfit)
                    << "    " << ConsensusFeaturePrinter(*cmit) << std::endl;
              }
              // We repeat the first feature handle at the end of the list.
              // This way you can generate closed line drawings
              // See Gnuplot set datafile commentschars
              consensus_elements_file << "L " << FeatureHandlePrinter(
                *cmit->begin()) << "    " << ConsensusFeaturePrinter(*cmit)
                  << std::endl;
            }
            consensus_elements_file.close();
          }

          // ----------------------------------------------------------------------

          if ( !consensus_features.empty() )
          {
            std::ofstream consensus_features_file(consensus_features.c_str());
            if ( !consensus_features_file )
            {
              throw Exception::UnableToCreateFile(__FILE__, __LINE__,
                __PRETTY_FUNCTION__, consensus_features);
            }

            std::map<Size,Size> map_id_to_map_num;
            std::vector<Size> map_num_to_map_id;
            std::vector<FeatureHandle> feature_handles;
            FeatureHandle feature_handle_NaN;
            feature_handle_NaN.setRT(std::numeric_limits<
                FeatureHandle::CoordinateType>::quiet_NaN());
            feature_handle_NaN.setMZ(std::numeric_limits<
                FeatureHandle::CoordinateType>::quiet_NaN());
            feature_handle_NaN.setIntensity(std::numeric_limits<
                FeatureHandle::IntensityType>::quiet_NaN());
            // feature_handle_NaN.setCharge(std::numeric_limits<Int>::max());

            for ( ConsensusMap::FileDescriptions::const_iterator fdit =
                consensus_map.getFileDescriptions().begin(); fdit
                != consensus_map.getFileDescriptions().end(); ++fdit )
            {
              map_id_to_map_num[fdit->first] = map_num_to_map_id.size();
              map_num_to_map_id.push_back(fdit->first);
            }

            consensus_features_file << "#  Consensus features extracted from "
                << in << " on " << date_time_now << std::endl
                << "# RT_cf MZ_cf Intensity_cf Charge_cf";
            for ( Size fhindex = 0; fhindex < map_num_to_map_id.size(); ++fhindex )
            {
              Size map_id = map_num_to_map_id[fhindex];
              consensus_features_file << "    RT_" << map_id << " MZ_"
                  << map_id << " Intensity_" << map_id << " Charge_" << map_id;
            }
            consensus_features_file << std::endl;

            for ( ConsensusMap::const_iterator cmit = consensus_map.begin(); cmit
                != consensus_map.end(); ++cmit )
            {
              {
                // please can anyone explain to me why putting the next two things into one statement doesnt work?
                std::vector<FeatureHandle> tmp(map_num_to_map_id.size(),
                  feature_handle_NaN);
                feature_handles.swap(tmp);
              }
              consensus_features_file << ConsensusFeaturePrinter(*cmit);
              for ( ConsensusFeature::const_iterator cfit = cmit->begin(); cfit
                  != cmit->end(); ++cfit )
              {
                feature_handles[map_id_to_map_num[cfit->getMapIndex()]] = *cfit;
              }
              for ( Size fhindex = 0; fhindex < feature_handles.size(); ++fhindex )
              {
                consensus_features_file << "    " << FeatureHandlePrinter(
                  feature_handles[fhindex]);
              }
              consensus_features_file << std::endl;
            }
            consensus_features_file.close();
          }

          // ----------------------------------------------------------------------

          if ( !out.empty() )
          {

            std::ofstream outstr(out.c_str());
            if ( !outstr )
            {
              throw Exception::UnableToCreateFile(__FILE__, __LINE__,
                __PRETTY_FUNCTION__, out);
            }

            outstr << "#  Consensus features extracted from " << in << " on "
                << date_time_now << std::endl;

            std::map<Size,Size> map_id_to_map_num;
            std::vector<Size> map_num_to_map_id;
            std::vector<FeatureHandle> feature_handles;
            FeatureHandle feature_handle_NaN;
            feature_handle_NaN.setRT(std::numeric_limits<
                FeatureHandle::CoordinateType>::quiet_NaN());
            feature_handle_NaN.setMZ(std::numeric_limits<
                FeatureHandle::CoordinateType>::quiet_NaN());
            feature_handle_NaN.setIntensity(std::numeric_limits<
                FeatureHandle::IntensityType>::quiet_NaN());
            feature_handle_NaN.setCharge(0); // just to be sure...
            // feature_handle_NaN.setCharge(std::numeric_limits<Int>::max()); // alternative ??

            // Its hard to predict which meta keys will be used in file descriptions.
            // So we assemble a list each time.  Represent keys by String, not UInt, for implicit sorting.
            std::set<String> all_file_desc_meta_keys;
            std::vector<UInt> tmp_meta_keys;
            for ( ConsensusMap::FileDescriptions::const_iterator fdit =
                consensus_map.getFileDescriptions().begin(); fdit
                != consensus_map.getFileDescriptions().end(); ++fdit )
            {
              map_id_to_map_num[fdit->first] = map_num_to_map_id.size();
              map_num_to_map_id.push_back(fdit->first);
              fdit->second.getKeys(tmp_meta_keys);
              for ( std::vector<UInt>::const_iterator kit =
                  tmp_meta_keys.begin(); kit != tmp_meta_keys.end(); ++kit )
              {
                all_file_desc_meta_keys.insert(
                  MetaInfoInterface::metaRegistry().getName(*kit));
              }
            }

            outstr << "#MAP" << sep << "id" << sep << "filename" << sep
                << "label" << sep << "size";
            for ( std::set<String>::const_iterator kit =
                all_file_desc_meta_keys.begin(); kit
                != all_file_desc_meta_keys.end(); ++kit )
            {
              outstr << sep << *kit;
            }
            outstr << std::endl;

            // list of maps (intentionally at the beginning, contrary to order in consensusXML)
            for ( ConsensusMap::FileDescriptions::const_iterator fdit =
                consensus_map.getFileDescriptions().begin(); fdit
                != consensus_map.getFileDescriptions().end(); ++fdit )
            {
              if ( no_ids )
              {
                outstr << "#";
              }
              outstr << "MAP" << sep << fdit->first << sep
                  << ( fdit->second.filename.empty() ? "\"\""
                      : fdit->second.filename ) << sep
                  << ( fdit->second.label.empty() ? "\"\"" : fdit->second.label )
                  << sep << fdit->second.size;
              for ( std::set<String>::const_iterator kit =
                  all_file_desc_meta_keys.begin(); kit
                  != all_file_desc_meta_keys.end(); ++kit )
              {
                if ( fdit->second.metaValueExists(*kit) )
                {
                  outstr << sep << fdit->second.getMetaValue(*kit);
                }
                else
                {
                  outstr << sep << "\"\""; // default is "" - reasonable?
                }
              }
              outstr << std::endl;
            }

            // stores one consensus feature per line
            if ( no_ids )
            {
              outstr << "#rt_cf" << sep << "mz_cf" << sep << "intensity_cf"
                  << sep << "charge_cf" << sep << "quality_cf";
              for ( Size fhindex = 0; fhindex < map_num_to_map_id.size(); ++fhindex )
              {
                Size map_id = map_num_to_map_id[fhindex];
                outstr << sep << "rt_" << map_id << sep << "mz_" << map_id
                    << sep << "intensity_" << map_id << sep << "charge_"
                    << map_id;
              }
              outstr << std::endl;
            }
            else
            {

              outstr << "#CONSENSUS" << sep << "rt_cf" << sep << "mz_cf" << sep
                  << "intensity_cf" << sep << "charge_cf" << sep << "quality_cf";
              for ( Size fhindex = 0; fhindex < map_num_to_map_id.size(); ++fhindex )
              {
                Size map_id = map_num_to_map_id[fhindex];
                outstr << sep << "rt_" << map_id << sep << "mz_" << map_id
                    << sep << "intensity_" << map_id << sep << "charge_"
                    << map_id;
              }
              outstr << std::endl;

              outstr << "#RUN" << sep << "RunID" << sep << "ScoreType" << sep
                  << "ScoreDirection" << sep << "Date/Time" << sep
                  << "SearchEngineVersion" << sep << "Parameters" << std::endl;

              outstr << "#PROTEIN" << sep << "Score" << sep << "Rank" << sep
                  << "Accession" << sep << "Sequence" << std::endl;

              outstr << "#UNASSIGNEDPEPTIDE" << sep << "rt" << sep << "mz"
                  << sep << "score" << sep << "rank" << sep << "sequence"
                  << sep << "charge" << sep << "AA_before" << sep << "AA_after"
                  << sep << "score_type" << sep << "search_identifier" <<
              // not sure if accessions really make sense, but it's copy&paste code anyway and predicted_RT is sensible
                  sep << "accessions" << sep << "predicted_RT" << std::endl;

              // protein accessions, predicted_RT

              outstr << "#PEPTIDE" << sep << "rt" << sep << "mz" << sep
                  << "score" << sep << "rank" << sep << "sequence" << sep
                  << "charge" << sep << "AA_before" << sep << "AA_after" << sep
                  << "score_type" << sep << "search_identifier" << std::endl;

            }

            //  proteins and unassigned peptides
            if ( !no_ids )
            {
              // proteins
              {
                for ( vector<ProteinIdentification>::const_iterator it =
                    consensus_map.getProteinIdentifications().begin(); it
                    != consensus_map.getProteinIdentifications().end(); ++it )
                {
                  String actual_id = it->getIdentifier();
                  // protein id header
                  outstr << "RUN" << sep << actual_id << sep
                      << it->getScoreType() << sep;
                  if ( it->isHigherScoreBetter() )
                  {
                    outstr << "higher-score-better" << sep;
                  }
                  else
                  {
                    outstr << "lower-score-better" << sep;
                  }
                  // using ISODate ensures that TOPP tests will run through regardless of locale setting
                  outstr
                      << it->getDateTime().toString(Qt::ISODate).toStdString()
                      << sep << it->getSearchEngineVersion() << sep;

                  // search parameters
                  ProteinIdentification::SearchParameters sp =
                      it->getSearchParameters();
                  outstr << "db=" << sp.db << ", db_version=" << sp.db_version
                      << ", taxonomy=" << sp.taxonomy << ", charges="
                      << sp.charges << ", mass_type=";
                  if ( sp.mass_type == ProteinIdentification::MONOISOTOPIC )
                  {
                    outstr << "monoisotopic";
                  }
                  else
                  {
                    outstr << "average";
                  }
                  outstr << ", fixed_modifications=";
                  for ( vector<String>::const_iterator mit =
                      sp.fixed_modifications.begin(); mit
                      != sp.fixed_modifications.end(); ++mit )
                  {
                    if ( mit != sp.fixed_modifications.begin() )
                    {
                      outstr << ";";
                    }
                    outstr << *mit;
                  }
                  outstr << ", variable_modifications=";
                  for ( vector<String>::const_iterator mit =
                      sp.variable_modifications.begin(); mit
                      != sp.variable_modifications.end(); ++mit )
                  {
                    if ( mit != sp.variable_modifications.begin() )
                    {
                      outstr << ";";
                    }
                    outstr << *mit;
                  }
                  outstr << ", enzyme=";
                  switch ( sp.enzyme )
                  {
                    case ProteinIdentification::TRYPSIN:
                      outstr << "Trypsin";
                      break;
                    case ProteinIdentification::PEPSIN_A:
                      outstr << "PepsinA";
                      break;
                    case ProteinIdentification::PROTEASE_K:
                      outstr << "ProteaseK";
                      break;
                    case ProteinIdentification::CHYMOTRYPSIN:
                      outstr << "ChymoTrypsin";
                      break;
                    default:
                      outstr << "unknown";
                  }
                  outstr << ", missed_cleavages=" << sp.missed_cleavages
                      << ", peak_mass_tolerance=" << sp.peak_mass_tolerance
                      << ", precursor_mass_tolerance="
                      << sp.precursor_tolerance << endl;

                  for ( vector<ProteinHit>::const_iterator pit =
                      it->getHits().begin(); pit != it->getHits().end(); ++pit )
                  {
                    outstr << "PROTEIN" << sep << pit->getScore() << sep
                        << pit->getRank() << sep << pit->getAccession() << sep
                        << pit->getSequence() << endl;
                  }
                }
              }

              // unassigned peptides
              {
                for ( vector<PeptideIdentification>::const_iterator pit =
                    consensus_map.getUnassignedPeptideIdentifications().begin(); pit
                    != consensus_map.getUnassignedPeptideIdentifications().end(); ++pit )
                {
                  for ( vector<PeptideHit>::const_iterator ppit =
                      pit->getHits().begin(); ppit != pit->getHits().end(); ++ppit )
                  {
                    outstr << "UNASSIGNEDPEPTIDE" << sep;
                    if ( pit->metaValueExists("RT") )
                    {
                      outstr << (DoubleReal) pit->getMetaValue("RT") << sep;
                    }
                    else
                    {
                      outstr << "-1" << sep;
                    }

                    if ( pit->metaValueExists("MZ") )
                    {
                      outstr << (DoubleReal) pit->getMetaValue("MZ") << sep;
                    }
                    else
                    {
                      outstr << "-1" << sep;
                    }
                    outstr << ppit->getScore() << sep << ppit->getRank() << sep
                        << ppit->getSequence() << sep << ppit->getCharge()
                        << sep << ppit->getAABefore() << sep
                        << ppit->getAAAfter() << sep << pit->getScoreType()
                        << sep << pit->getIdentifier() << sep;

                    for ( vector<String>::const_iterator ait =
                        ppit->getProteinAccessions().begin(); ait
                        != ppit->getProteinAccessions().end(); ++ait )
                    {
                      if ( ait != ppit->getProteinAccessions().begin() )
                      {
                        outstr << ";";
                      }
                      outstr << *ait;
                    }
                    if ( ppit->metaValueExists("predicted_RT") )
                    {
                      outstr << sep << ppit->getMetaValue("predicted_RT");
                    }
                    else
                    {
                      outstr << sep << "-1";
                    }
                    /* first_dim_... stuff not supported for now */
                    //  if ( first_dim_rt )
                    //  {
                    //    if ( pit->metaValueExists("first_dim_rt") )
                    //    {
                    //      outstr << sep << pit->getMetaValue("first_dim_rt");
                    //    }
                    //    else
                    //    {
                    //      outstr << sep << "-1";
                    //    }
                    //    if ( ppit->metaValueExists("predicted_RT_first_dim") )
                    //    {
                    //      outstr << sep << ppit->getMetaValue(
                    //        "predicted_RT_first_dim");
                    //    }
                    //    else
                    //    {
                    //      outstr << sep << "-1";
                    //    }
                    //  }
                    outstr << std::endl;
                  }
                }
              }
            }

            for ( ConsensusMap::const_iterator cmit = consensus_map.begin(); cmit
                != consensus_map.end(); ++cmit )
            {
              {
                // please can anyone explain to me why putting the next two things into one statement doesnt work?
                std::vector<FeatureHandle> tmp(map_num_to_map_id.size(),
                  feature_handle_NaN);
                feature_handles.swap(tmp);
              }
              if ( !no_ids ) outstr << "CONSENSUS" << sep;
              outstr << precisionWrapper(cmit->getRT()) << sep
                  << precisionWrapper(cmit->getMZ()) << sep
                  << precisionWrapper(cmit->getIntensity()) << sep
                  << cmit->getCharge() << sep << cmit->getQuality();
              for ( ConsensusFeature::const_iterator cfit = cmit->begin(); cfit
                  != cmit->end(); ++cfit )
              {
                feature_handles[map_id_to_map_num[cfit->getMapIndex()]] = *cfit;
              }
              for ( Size fhindex = 0; fhindex < feature_handles.size(); ++fhindex )
              {
                outstr << sep << precisionWrapper(
                  feature_handles[fhindex].getRT()) << sep << precisionWrapper(
                  feature_handles[fhindex].getMZ()) << sep << precisionWrapper(
                  feature_handles[fhindex].getIntensity()) << sep
                    << feature_handles[fhindex].getCharge();
              }
              outstr << std::endl;

              // peptide ids
              if ( !no_ids )
              {
                for ( vector<PeptideIdentification>::const_iterator pit =
                    cmit->getPeptideIdentifications().begin(); pit
                    != cmit->getPeptideIdentifications().end(); ++pit )
                {
                  for ( vector<PeptideHit>::const_iterator ppit =
                      pit->getHits().begin(); ppit != pit->getHits().end(); ++ppit )
                  {
                    outstr << "PEPTIDE" << sep;
                    if ( pit->metaValueExists("RT") )
                    {
                      outstr << precisionWrapper(pit->getMetaValue("RT"))
                          << sep;
                    }
                    else
                    {
                      outstr << "-1" << sep;
                    }

                    if ( pit->metaValueExists("MZ") )
                    {
                      outstr << precisionWrapper(pit->getMetaValue("MZ"))
                          << sep;
                    }
                    else
                    {
                      outstr << "-1" << sep;
                    }
                    outstr << precisionWrapper(ppit->getScore()) << sep
                        << precisionWrapper(ppit->getRank()) << sep
                        << ppit->getSequence() << sep << ppit->getCharge()
                        << sep << ppit->getAABefore() << sep
                        << ppit->getAAAfter() << sep << pit->getScoreType()
                        << sep << pit->getIdentifier() << std::endl;
                  }
                }
              }
            }
          }

          return EXECUTION_OK;
        }
        else if ( in_type == FileTypes::IDXML )
        {
          vector<ProteinIdentification> prot_ids;
          vector<PeptideIdentification> pep_ids;
          String document_id;
          IdXMLFile().load(in, prot_ids, pep_ids, document_id);

          counter = 0;
          ofstream txt_out(out.c_str());

          txt_out << "#RUN" << sep << "RunID" << sep << "ScoreType" << sep
              << "ScoreDirection" << sep << "Date/Time" << sep
              << "SearchEngineVersion" << sep << "Parameters" << endl;
          txt_out << "#PROTEIN" << sep << "Score" << sep << "Rank" << sep
              << "Accession" << sep << "Sequence" << endl;

          if ( first_dim_rt )
          {
            txt_out << "#PEPTIDE" << sep << "RT" << sep << "MZ" << sep
                << "Score" << sep << "Rank" << sep << "Sequence" << sep
                << "Charge" << sep << "AABefore" << sep << "AAAfter" << sep
                << "Accessions" << sep << "predicted_RT" << sep
                << "RT_first_dim" << sep << "predicted_RT_first_dim" << endl;
          }
          else
          {
            txt_out << "#PEPTIDE" << sep << "RT" << sep << "MZ" << sep
                << "Score" << sep << "Rank" << sep << "Sequence" << sep
                << "Charge" << sep << "AABefore" << sep << "AAAfter" << sep
                << "Accessions" << sep << "predicted_RT" << endl;
          }

          for ( vector<ProteinIdentification>::const_iterator it =
              prot_ids.begin(); it != prot_ids.end(); ++it )
          {
            String actual_id = it->getIdentifier();
            if ( !getFlag_("peptides_only") )
            {
              // protein id header
              txt_out << "RUN" << sep << actual_id << sep << it->getScoreType()
                  << sep;
              if ( it->isHigherScoreBetter() )
              {
                txt_out << "higher-score-better" << sep;
              }
              else
              {
                txt_out << "lower-score-better" << sep;
              }
              // using ISODate ensures that TOPP tests will run through regardless of locale setting
              txt_out << it->getDateTime().toString(Qt::ISODate).toStdString()
                  << sep << it->getSearchEngineVersion() << sep;

              // search parameters
              ProteinIdentification::SearchParameters sp =
                  it->getSearchParameters();
              txt_out << "db=" << sp.db << ", db_version=" << sp.db_version
                  << ", taxonomy=" << sp.taxonomy << ", charges=" << sp.charges
                  << ", mass_type=";
              if ( sp.mass_type == ProteinIdentification::MONOISOTOPIC )
              {
                txt_out << "monoisotopic";
              }
              else
              {
                txt_out << "average";
              }
              txt_out << ", fixed_modifications=";
              for ( vector<String>::const_iterator mit =
                  sp.fixed_modifications.begin(); mit
                  != sp.fixed_modifications.end(); ++mit )
              {
                if ( mit != sp.fixed_modifications.begin() )
                {
                  txt_out << ";";
                }
                txt_out << *mit;
              }
              txt_out << ", variable_modifications=";
              for ( vector<String>::const_iterator mit =
                  sp.variable_modifications.begin(); mit
                  != sp.variable_modifications.end(); ++mit )
              {
                if ( mit != sp.variable_modifications.begin() )
                {
                  txt_out << ";";
                }
                txt_out << *mit;
              }
              txt_out << ", enzyme=";
              switch ( sp.enzyme )
              {
                case ProteinIdentification::TRYPSIN:
                  txt_out << "Trypsin";
                  break;
                case ProteinIdentification::PEPSIN_A:
                  txt_out << "PepsinA";
                  break;
                case ProteinIdentification::PROTEASE_K:
                  txt_out << "ProteaseK";
                  break;
                case ProteinIdentification::CHYMOTRYPSIN:
                  txt_out << "ChymoTrypsin";
                  break;
                default:
                  txt_out << "unknown";
              }
              txt_out << ", missed_cleavages=" << sp.missed_cleavages
                  << ", peak_mass_tolerance=" << sp.peak_mass_tolerance
                  << ", precursor_mass_tolerance=" << sp.precursor_tolerance
                  << endl;

              for ( vector<ProteinHit>::const_iterator pit =
                  it->getHits().begin(); pit != it->getHits().end(); ++pit )
              {
                txt_out << "PROTEIN" << sep << pit->getScore() << sep
                    << pit->getRank() << sep << pit->getAccession() << sep
                    << pit->getSequence() << endl;
              }
            }

            if ( !getFlag_("proteins_only") )
            {
              // slight improvement on big idXML files with many different runs:
              // index the identifiers and peptide ids to avoid running over them
              // again and again
              for ( vector<PeptideIdentification>::const_iterator pit =
                  pep_ids.begin(); pit != pep_ids.end(); ++pit )
              {
                if ( pit->getIdentifier() == actual_id )
                {
                  for ( vector<PeptideHit>::const_iterator ppit =
                      pit->getHits().begin(); ppit != pit->getHits().end(); ++ppit )
                  {
                    txt_out << "PEPTIDE" << sep;
                    if ( pit->metaValueExists("RT") )
                    {
                      txt_out << (DoubleReal) pit->getMetaValue("RT") << sep;
                    }
                    else
                    {
                      txt_out << "-1" << sep;
                    }

                    if ( pit->metaValueExists("MZ") )
                    {
                      txt_out << (DoubleReal) pit->getMetaValue("MZ") << sep;
                    }
                    else
                    {
                      txt_out << "-1" << sep;
                    }
                    txt_out << ppit->getScore() << sep << ppit->getRank()
                        << sep << ppit->getSequence() << sep
                        << ppit->getCharge() << sep << ppit->getAABefore()
                        << sep << ppit->getAAAfter() << sep;

                    for ( vector<String>::const_iterator ait =
                        ppit->getProteinAccessions().begin(); ait
                        != ppit->getProteinAccessions().end(); ++ait )
                    {
                      if ( ait != ppit->getProteinAccessions().begin() )
                      {
                        txt_out << ";";
                      }
                      txt_out << *ait;
                    }
                    if ( ppit->metaValueExists("predicted_RT") )
                    {
                      txt_out << sep << ppit->getMetaValue("predicted_RT");
                    }
                    else
                    {
                      txt_out << sep << "-1";
                    }
                    if ( first_dim_rt )
                    {
                      if ( pit->metaValueExists("first_dim_rt") )
                      {
                        txt_out << sep << pit->getMetaValue("first_dim_rt");
                      }
                      else
                      {
                        txt_out << sep << "-1";
                      }
                      if ( ppit->metaValueExists("predicted_RT_first_dim") )
                      {
                        txt_out << sep << ppit->getMetaValue(
                          "predicted_RT_first_dim");
                      }
                      else
                      {
                        txt_out << sep << "-1";
                      }
                    }
                    txt_out << endl;
                  }
                }
              }
            }
          }

          txt_out.close();
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
