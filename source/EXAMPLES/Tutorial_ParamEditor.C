#include <QtGui/QApplication>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

using namespace OpenMS;
using namespace std;

Int main(int argc, const char** argv)
{	
	QApplication app(argc,const_cast<char**>(argv));
	
	Param param;
	param.load(argv[1]);
	
	ParamEditor* editor = new ParamEditor(0);
	editor->loadEditable(param);
	editor->show();
	
	app.exec();
	
	editor->store();
	param.store(argv[1]);
	
	return 0;
} //end of main
