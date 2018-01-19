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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <fstream>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_OpenSwathRewriteToFeatureXML OpenSwathRewriteToFeatureXML

  @brief Combines featureXML and mProphet tsv to FDR filtered featureXML.

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_OpenSwathRewriteToFeatureXML.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_OpenSwathRewriteToFeatureXML.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPOpenSwathRewriteToFeatureXML : 
  public TOPPBase, 
  public ProgressLogger
{
 public:

  TOPPOpenSwathRewriteToFeatureXML()
    : TOPPBase("OpenSwathRewriteToFeatureXML","Combines featureXML and mProphet tsv to FDR filtered featureXML.", false)
  {
  }

 protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("csv","<file>","","mProphet tsv output file: \"all_peakgroups.xls\"", false);
    setValidFormats_("csv", ListUtils::create<String>("csv"));
    
    registerInputFile_("featureXML","<file>","","input featureXML file");
    setValidFormats_("featureXML", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out","<file>","","output featureXML file");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));

    registerDoubleOption_("FDR_cutoff", "<double>", -1, "FDR cutoff (e.g. to remove all features with a an m_score above 0.05 use 0.05 here)", false);
  }

  void applyFDRcutoff(FeatureMap & feature_map, double cutoff, String fdr_name)
  {
    FeatureMap out_feature_map = feature_map;
    out_feature_map.clear(false);
    for (Size i = 0; i < feature_map.size(); i++)
    {
      if ((double)feature_map[i].getMetaValue(fdr_name) < cutoff)
      {
        out_feature_map.push_back(feature_map[i]);
      }
    }
    feature_map = out_feature_map;
  }

  void processInput(const char * filename, FeatureMap & feature_map)
  {
    FeatureMap out_feature_map = feature_map;
    std::map<String, int> added_already;
    out_feature_map.clear(false);

    std::map<String, Feature*> feature_map_ref;
    //for (FeatureMap::iterator feature = feature_map.begin(); feature != feature_map.end(); feature++)
    for (Size i = 0; i < feature_map.size(); i++)
    {
      feature_map_ref[feature_map[i].getUniqueId()] = &feature_map[i];
    }

    std::ifstream data(filename);
    std::string   line;

    // Read header
    std::getline(data, line);
    // std::map<int, String> header_dict; // not used
    std::map<String, int> header_dict_inv;
    {
      std::stringstream          lineStream(line);
      std::string                cell;
      int cnt = 0;
      while (std::getline(lineStream,cell,'\t'))
      {
        //header_dict[cnt] = cell;
        header_dict_inv[cell] = cnt;
        cnt++;
      }
    }

    if (header_dict_inv.find("id") == header_dict_inv.end() || 
        header_dict_inv.find("m_score") == header_dict_inv.end() || 
        header_dict_inv.find("d_score") == header_dict_inv.end() )
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: The tsv file is expected to have at least the following headers: id, m_score, d_score. " );
    }

    // Read file
    std::vector<std::string> current_row;
    std::string                cell;
    int line_nr = 0;
    double m_score, d_score;
    while (std::getline(data, line))
    {
      line_nr++;
      current_row.clear();
      std::stringstream  lineStream(line);
      while (std::getline(lineStream,cell,'\t'))
      {
        current_row.push_back(cell);
      }

      String id = current_row[header_dict_inv["id"]];
      id = id.substitute("f_", ""); 
      try
      {
        m_score = ((String)current_row[header_dict_inv["m_score"]]).toDouble();
      }
      catch (char* /*str*/)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: Could not convert String" + ((String)current_row[header_dict_inv["m_score"]]) + " on line " + String(line_nr));
      }
      try
      {
        d_score = ((String)current_row[header_dict_inv["d_score"]]).toDouble();
      }
      catch (char* /*str*/)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: Could not convert String" + ((String)current_row[header_dict_inv["d_score"]]) + " on line " + String(line_nr));
      }

      if (feature_map_ref.find(id) != feature_map_ref.end() )
      {
        Feature* feature = feature_map_ref.find(id)->second;
        feature->setMetaValue("m_score", m_score);
        feature->setMetaValue("d_score", d_score);
        // we are not allowed to have duplicate unique ids
        if (added_already.find(id) != added_already.end())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: Duplicate id found in CSV file: " + id );
        }
        out_feature_map.push_back(*feature);
      }
    }
    feature_map = out_feature_map;
  }

  ExitCodes main_(int , const char**) override
  {

  String feature_file = getStringOption_("featureXML");
  String csv = getStringOption_("csv");
  String out = getStringOption_("out");
  double fdr_cutoff = getDoubleOption_("FDR_cutoff");

  FeatureMap feature_map;
  FeatureXMLFile().load(feature_file, feature_map);

  if (csv.size() > 0)
  {
    processInput(csv.c_str(), feature_map);
  }

  if (fdr_cutoff >= 0)
  {
    applyFDRcutoff(feature_map, fdr_cutoff, "m_score");
  }


  feature_map.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
  FeatureXMLFile().store(out, feature_map);

  return EXECUTION_OK;

  }

};

int main( int argc, const char** argv )
{

  TOPPOpenSwathRewriteToFeatureXML tool;
  int code = tool.main(argc,argv);
  return code;

}

/// @endcond

