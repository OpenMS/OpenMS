#include <QtGui/QApplication>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>

using namespace OpenMS;
using namespace std;

Int main(int argc, const char** argv)
{  
  QApplication app(argc,const_cast<char**>(argv));
  
  PeakMap exp;
  exp.resize(1);
  DTAFile().load("data/Tutorial_Spectrum1D.dta",exp[0]);
  
  Spectrum1DWidget* widget = new Spectrum1DWidget(Param(),0);
  widget->canvas()->addLayer(exp);
  widget->show();
  
  return app.exec();
} //end of main
