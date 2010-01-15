#include<QtGui/QApplication>
#include<OpenMS/VISUAL/ListEditor.h>

using namespace OpenMS;

int main(int argc, char** argv)
{
  QApplication app(argc, argv);
  ListEditor *listeditor = new ListEditor;
  listeditor->show();
  return app.exec();  
  
  return 0;
}
