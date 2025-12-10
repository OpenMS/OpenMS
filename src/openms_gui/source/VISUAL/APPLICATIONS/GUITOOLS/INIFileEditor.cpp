// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#include <OpenMS/VISUAL/APPLICATIONS/INIFileEditorWindow.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>

// Qt
#include <QSurfaceFormat>
#if !defined(__APPLE__)
  #include <QtWidgets/QStyleFactory>
#endif

#ifdef OPENMS_WINDOWSPLATFORM
#   ifndef _WIN32_WINNT
#       define _WIN32_WINNT 0x0501 // Win XP (and above)
#   endif
#   include <Windows.h>
#endif

using namespace OpenMS;
using namespace std;

/**
    @page TOPP_INIFileEditor INIFileEditor

    @brief Can be used to visually edit INI files of TOPP tools.

    The values can be edited by double-clicking or pressing F2.

    The documentation of each value is shown in the text area on the bottom of the widget.

    @image html INIFileEditor.png
*/

int main(int argc, const char** argv)
{
#ifdef OPENMS_WINDOWSPLATFORM
  qputenv("QT_QPA_PLATFORM", "windows:darkmode=0"); // disable dark mode on Windows, since our buttons etc are not designed for it
#endif

  std::map<std::string, std::string> options, flags, option_lists;
  options["-print"] = "print";
  flags["--help"] = "help";
  Param param;
  param.parseCommandLine(argc, argv, options, flags, option_lists);

  //catch command line errors
  if (param.exists("help") //help requested
     || argc > 3 //too many arguments
     || (argc == 3 && !param.exists("print")) //three argument but no -print
     || (param.exists("print") && param.getValue("print") == "") //-print but no file given
      )
  {
    cerr << endl
         << "INIFileEditor -- An editor for OpenMS configuration files." << endl
         << endl
         << "Usage:" << endl
         << " INIFileEditor [options] [file]" << endl
         << endl
         << "Options are:" << endl
         << " --help         Shows this help and exits" << endl
         << " -print <file>  Prints the content of the file to the command line and exits" << endl
         << endl;
    return 0;
  }

  //print a ini file as text
  if (param.exists("print"))
  {
    Param data;
    ParamXMLFile paramFile;
    try
    {
      paramFile.load(param.getValue("print").toString(), data);
      for (Param::ParamIterator it = data.begin(); it != data.end(); ++it)
      {
        cout << it.getName() << " = " << it->value << endl;
      }
    }
    catch (Exception::BaseException& e)
    {
      OPENMS_LOG_ERROR << "Error while parsing file '" << param.getValue("print") << "'\n";
      OPENMS_LOG_ERROR << e << "\n";
    }

    return 0;
  }

  //Create window
#if defined(__APPLE__)
  // see https://bugreports.qt.io/browse/QTBUG-104871
  // if you link to QtWebEngine and the corresponding macros are enabled, it will
  // try to default to OpenGL 4.1 on macOS (for hardware acceleration of WebGL in Chromium, which we do not need yet)
  // but our OpenGL code for 3D View is written in OpenGL 2.x.
  // Now we force 2.1 which is also available on all? Macs.
  QSurfaceFormat format;
  format.setVersion(2, 1); // the default is 2, 0
  QSurfaceFormat::setDefaultFormat(format); // should be done before creating a QApplication
#endif

  QApplicationTOPP app(argc, const_cast<char**>(argv));

  INIFileEditorWindow editor_window;

  //Open passed file
  if (argc == 2)
  {
    //cout << "OPEN: "  << argv[1] << endl;
    editor_window.openFile(argv[1]);
  }

#ifdef OPENMS_WINDOWSPLATFORM
  FreeConsole(); // get rid of console window at this point (we will not see any console output from this point on)
  AttachConsole(-1); // if the parent is a console, reattach to it - so we can see debug output - a normal user will usually not use cmd.exe to start a GUI)
#endif

  editor_window.show();
  return app.exec();
}
