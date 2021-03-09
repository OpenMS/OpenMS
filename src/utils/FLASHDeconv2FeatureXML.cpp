//
// Created by JihyungKim on 11/30/20.
//
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <fstream>
#include <OpenMS/CONCEPT/LogStream.h>
#include <boost/algorithm/string.hpp>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/Constants.h>

using namespace OpenMS;
using namespace std;

class TOPPFLASHDeconv2FeatureXML :
    public TOPPBase
{
public:
  TOPPFLASHDeconv2FeatureXML():
      TOPPBase("FLASHDeconv2FeatureXML", "Convert FLASHDeconv result file to FeatureXML format", false, {}, false)
  {
  }
protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file", true);
    setValidFormats_("in", ListUtils::create<String>("tsv"));
    registerIntOption_("num_feat", "<integer>", -1, "the number of features(default: -1, full features in input file", false);
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_path = getStringOption_("in");
    int number_of_features = getIntOption_("num_feat");
    if (number_of_features==-1)
    {
      int numLines = 0;
      std::ifstream in(inputfile_path);
      std::string unused;
      while ( std::getline(in, unused) )
        ++numLines;
      number_of_features = numLines-1; // exclude the header
    }

    // output file
    String outfilePath = inputfile_path.substr(0, inputfile_path.find_last_of(".")) + ".FD.featureXML";
    OPENMS_LOG_INFO << "output file : " << outfilePath << endl;
    FeatureMap fMap;
    fMap.reserve(number_of_features);

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    ifstream infile(inputfile_path);
    String line;

    // prepare header vector
    vector<String> header_vec;
    getline(infile, line); // header line
    boost::split(header_vec, line, boost::is_any_of("\t"));
    auto mass_column_index = distance(header_vec.begin(), find(header_vec.begin(), header_vec.end(), "MonoisotopicMass"));
    auto startRT_column_index = distance(header_vec.begin(), find(header_vec.begin(), header_vec.end(), "StartRetentionTime"));
    auto endRT_column_index = distance(header_vec.begin(), find(header_vec.begin(), header_vec.end(), "EndRetentionTime"));
    auto apexRT_column_index = distance(header_vec.begin(), find(header_vec.begin(), header_vec.end(), "ApexRetentionTime"));
    auto startCS_column_index = distance(header_vec.begin(), find(header_vec.begin(), header_vec.end(), "MinCharge"));
    auto endCS_column_index = distance(header_vec.begin(), find(header_vec.begin(), header_vec.end(), "MaxCharge"));
    auto sumInt_column_index = distance(header_vec.begin(), find(header_vec.begin(), header_vec.end(), "SumIntensity"));
    auto isoCoScore_column_index = distance(header_vec.begin(), find(header_vec.begin(), header_vec.end(), "IsotopeCosineScore"));
    auto perIsoInt_column_index = distance(header_vec.begin(), find(header_vec.begin(), header_vec.end(), "PerIsotopeIntensity"));

    // read file
    for (Size numOfF = 0; numOfF<number_of_features; numOfF++) {
      std::getline(infile, line);
      vector<String> results;
      boost::split(results, line, boost::is_any_of("\t"));

      auto isoStr = results[perIsoInt_column_index];
      auto mass = stod(results[mass_column_index]);
      auto startRT = stod(results[startRT_column_index]);
      auto endRT = stod(results[endRT_column_index]);
      auto apexRT = stod(results[apexRT_column_index]);
      auto sumInty = stod(results[sumInt_column_index]);
      auto numOfIso = count(isoStr.begin(), isoStr.end(), ';');
      auto score = stod(results[isoCoScore_column_index]);

      OPENMS_LOG_INFO << to_string(mass) << endl;

      auto cs = stoi(results[startCS_column_index]);
      auto cs_end = stoi(results[endCS_column_index]) +1;
      for (; cs<cs_end; cs++){
        Feature f;
        // set label
        f.setCharge(cs);
        f.setOverallQuality(score);
        f.setRT(apexRT);
        f.setIntensity(sumInty);

        // set apex mz
        auto monoMZ = (mass + cs*Constants::PROTON_MASS_U) / cs;
        f.setMZ(monoMZ);

        // get isotopes
        for (auto iso = 0; iso<numOfIso; iso++){
          auto isoMZ = Constants::C13C12_MASSDIFF_U / cs * iso + monoMZ;
          // draw convexHulls
          ConvexHull2D::PointArrayType hull_points(4);
          hull_points[0][0] = startRT; // min RT
          hull_points[0][1] = isoMZ; // MZ
          hull_points[1][0] = endRT; // max RT
          hull_points[1][1] = isoMZ; //MZ
          // redunant values for viewer
          hull_points[2][0] = startRT; // min RT
          hull_points[2][1] = isoMZ; // MZ
          hull_points[3][0] = endRT; // max RT
          hull_points[3][1] = isoMZ; //MZ

          ConvexHull2D hull;
          hull.addPoints(hull_points);

          f.getConvexHulls().push_back(hull);
        }
        fMap.push_back(f);
      }

    }


    FeatureXMLFile().store(outfilePath, fMap);

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFLASHDeconv2FeatureXML tool;
  return tool.main(argc, argv);
}
