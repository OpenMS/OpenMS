//
// Created by JihyungKim on 1/29/21.
//
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

class TOPPFeatureFinderIntact2GraphML :
    public TOPPBase
{
public :
  TOPPFeatureFinderIntact2GraphML():
    TOPPBase("FeatureFinderIntact2GraphML", "", false, {}, false)
  {}

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file", true);
    setValidFormats_("in", ListUtils::create<String>("tsv"));
  }

};