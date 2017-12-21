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
// $Authors: Marc Sturm, Clemens Groepl, Steffen Sass $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>

#include <OpenMS/KERNEL/ConversionHelper.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_FeatureLinkerBase FeatureLinkerBase

    @brief Base class for different FeatureLinker tools.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureLinkerBase :
  public TOPPBase, public ProgressLogger
{

public:
  TOPPFeatureLinkerBase(String name, String description, bool official = true) :
    TOPPBase(name, description, official)
  {
  }

protected:
  void registerOptionsAndFlags_() override   // only for "unlabeled" algorithms!
  {
    registerInputFileList_("in", "<files>", ListUtils::create<String>(""), "input files separated by blanks", true);
    setValidFormats_("in", ListUtils::create<String>("featureXML,consensusXML"));
    registerOutputFile_("out", "<file>", "", "Output file", true);
    setValidFormats_("out", ListUtils::create<String>("consensusXML"));
    registerInputFile_("design", "<file>", "", "input file containing the experimental design", false);
    setValidFormats_("design", ListUtils::create<String>("tsv"));
    addEmptyLine_();
    registerFlag_("keep_subelements", "For consensusXML input only: If set, the sub-features of the inputs are transferred to the output.");
  }

  ExitCodes common_main_(FeatureGroupingAlgorithm * algorithm,
                         bool labeled = false)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    StringList ins;
    if (labeled)
    {
      ins.push_back(getStringOption_("in"));
    }
    else
    {
      ins = getStringList_("in");
    }
    String out = getStringOption_("out");
    
    //-------------------------------------------------------------
    // check for valid input
    //-------------------------------------------------------------
    // check if all input files have the correct type
    FileTypes::Type file_type = FileHandler::getType(ins[0]);
    for (Size i = 0; i < ins.size(); ++i)
    {
      if (FileHandler::getType(ins[i]) != file_type)
      {
        writeLog_("Error: All input files must be of the same type!");
        return ILLEGAL_PARAMETERS;
      }
    }

    //-------------------------------------------------------------
    // set up algorithm
    //-------------------------------------------------------------
    Param algorithm_param = getParam_().copy("algorithm:", true);
    writeDebug_("Used algorithm parameters", algorithm_param, 3);
    algorithm->setParameters(algorithm_param);

    //-------------------------------------------------------------
    // perform grouping
    //-------------------------------------------------------------
    // load input
    ConsensusMap out_map;
    StringList ms_run_locations;

    String design_file;

    // TODO: support design in labeled feature linker
    if (!labeled)
    {
      design_file = getStringOption_("design");
    }

    if (file_type == FileTypes::CONSENSUSXML && !design_file.empty())
    {
      writeLog_("Error: Using fractionated design with consensusXML als input is not supported!");
      return ILLEGAL_PARAMETERS;
    }
  
  
    if (file_type == FileTypes::FEATUREXML)
    {
      //-------------------------------------------------------------
      // Extract (optional) fraction identifiers and associate with featureXMLs
      //-------------------------------------------------------------

      // determine map of fractions to runs
      map<unsigned, set<unsigned> > frac2run;

      if (!design_file.empty())
      {
        // parse design file and determine fractions
        ExperimentalDesign ed;
        ExperimentalDesign().load(design_file, ed);

        // determine if design defines more than one fraction
        frac2run = ed.getFractionToRunsMapping();

        // check if all fractions have the same number of MS runs associated
        if (!ed.sameNrOfRunsPerFraction())
        {
          writeLog_("Error: Number of runs must match for every fraction!");
          return ILLEGAL_PARAMETERS;
        }
      }
      else // no design file given
      {
        for (Size i = 0; i != ins.size(); ++i)
        {
          frac2run[1].insert(i + 1); // associate each run with fraction 1
        }
      }

      vector<ConsensusMap > maps(ins.size());
      FeatureXMLFile f;
      FeatureFileOptions param = f.getOptions();
      // to save memory don't load convex hulls and subordinates
      param.setLoadSubordinates(false);
      param.setLoadConvexHull(false);
      f.setOptions(param);

      Size progress = 0;
      setLogType(ProgressLogger::CMD);
      startProgress(0, ins.size(), "reading input");
      for (Size i = 0; i < ins.size(); ++i)
      {
        FeatureMap tmp;
        f.load(ins[i], tmp);
        out_map.getFileDescriptions()[i].filename = ins[i];
        out_map.getFileDescriptions()[i].size = tmp.size();
        out_map.getFileDescriptions()[i].unique_id = tmp.getUniqueId();

        // copy over information on the primary MS run
        StringList ms_runs;
        tmp.getPrimaryMSRunPath(ms_runs);
        ms_run_locations.insert(ms_run_locations.end(), ms_runs.begin(), ms_runs.end());

        // to save memory, remove convex hulls, subordinates:
        for (FeatureMap::Iterator it = tmp.begin(); it != tmp.end();
             ++it)
        {
          it->getSubordinates().clear();
          it->getConvexHulls().clear();
          it->clearMetaInfo();
        }

        MapConversion::convert(i, tmp, maps[i]);

        maps[i].updateRanges();

        setProgress(progress++);
      }
      endProgress();

      // exception for "labeled" algorithms: copy file descriptions
      if (labeled)
      {
        out_map.getFileDescriptions()[1] = out_map.getFileDescriptions()[0];
        out_map.getFileDescriptions()[0].label = "light";
        out_map.getFileDescriptions()[1].label = "heavy";
      }

      ////////////////////////////////////////////////////
      // invoke feature grouping algorithm
      
      if (frac2run.size() == 1) // group one fraction
      {
        algorithm->group(maps, out_map);
      }
      else // group multiple fractions
      {
        writeDebug_(String("Grouping ") + String(frac2run.size()) + " fractions.", 3);
        writeDebug_(String("Stored in ") + String(maps.size()) + " maps.", 3);
        for (Size i = 1; i <= frac2run.size(); ++i)
        {
          vector<ConsensusMap> fraction_maps;
          for (set<unsigned>::const_iterator sit = frac2run[i].begin(); sit != frac2run[i].end(); ++sit)
          {
            fraction_maps.push_back(maps[*sit - 1]); // TODO: *sit is currently the run identifier but we need to know the corresponding feature index in ins
          }
          algorithm->group(fraction_maps, out_map);
        }
      }
    }
    else
    {
      vector<ConsensusMap> maps(ins.size());
      ConsensusXMLFile f;
      for (Size i = 0; i < ins.size(); ++i)
      {
        f.load(ins[i], maps[i]);
        maps[i].updateRanges();
        // copy over information on the primary MS run
        StringList ms_runs;
        maps[i].getPrimaryMSRunPath(ms_runs);
        ms_run_locations.insert(ms_run_locations.end(), ms_runs.begin(), ms_runs.end());
      }
      // group
      algorithm->group(maps, out_map);

      // set file descriptions:
      bool keep_subelements = getFlag_("keep_subelements");
      if (!keep_subelements)
      {
        for (Size i = 0; i < ins.size(); ++i)
        {
          out_map.getFileDescriptions()[i].filename = ins[i];
          out_map.getFileDescriptions()[i].size = maps[i].size();
          out_map.getFileDescriptions()[i].unique_id = maps[i].getUniqueId();
        }
      }
      else
      {
        // components of the output map are not the input maps themselves, but
        // the components of the input maps:
        algorithm->transferSubelements(maps, out_map);
      }
    }

    // assign unique ids
    out_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    // annotate output with data processing info
    addDataProcessing_(out_map,
                       getProcessingInfo_(DataProcessing::FEATURE_GROUPING));

    // set primary MS runs
    out_map.setPrimaryMSRunPath(ms_run_locations);

    // sort list of peptide identifications in each consensus feature by map index
    out_map.sortPeptideIdentificationsByMapIndex();

    // write output
    ConsensusXMLFile().store(out, out_map);

    // some statistics
    map<Size, UInt> num_consfeat_of_size;
    for (ConsensusMap::const_iterator cmit = out_map.begin();
         cmit != out_map.end(); ++cmit)
    {
      ++num_consfeat_of_size[cmit->size()];
    }

    LOG_INFO << "Number of consensus features:" << endl;
    for (map<Size, UInt>::reverse_iterator i = num_consfeat_of_size.rbegin();
         i != num_consfeat_of_size.rend(); ++i)
    {
      LOG_INFO << "  of size " << setw(2) << i->first << ": " << setw(6) 
               << i->second << endl;
    }
    LOG_INFO << "  total:      " << setw(6) << out_map.size() << endl;

    return EXECUTION_OK;
  }

};

/// @endcond
