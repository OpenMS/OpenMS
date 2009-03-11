#include <QtGui/QApplication>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

using namespace OpenMS;
using namespace std;

Int main(int argc, const char** argv)
{  
  QApplication app(argc,const_cast<char**>(argv));
  
  Param param;
  param.load("data/Tutorial_ParamEditor.ini");
  
  ParamEditor* editor = new ParamEditor(0);
  editor->load(param);
  editor->show();
  
  app.exec();
  
  editor->store();
  param.store("output/Tutorial_ParamEditor_out.ini");
  
  return 0;
} //end of main
