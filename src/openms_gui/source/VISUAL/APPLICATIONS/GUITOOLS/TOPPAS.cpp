// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

/**
    @page TOPP_TOPPAS TOPPAS

    @brief An assistant for GUI-driven TOPP workflow design.

  TOPPAS allows to create, edit, open, save, and run TOPP workflows. Pipelines
  can be created conveniently in a GUI by means of mouse interactions. The
  parameters of all involved tools can be edited within the application
  and are also saved as part of the pipeline definition in the @em .toppas file.
  Furthermore, TOPPAS interactively performs validity checks during the pipeline
  editing process, in order to make it more difficult to create an invalid workflow.
  Once set up and saved, a workflow can also be run without the GUI using
  the @em ExecutePipeline TOPP tool.

  The following figure shows a simple example pipeline that has just been created
  and executed successfully:

  @image html TOPPAS_simple_example.png

  More information about TOPPAS can be found in the @ref TOPPAS_tutorial.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_TOPPAS.cli
*/

//QT
#include <QApplication>
#include <QPainter>
#include <QtWidgets/QSplashScreen>
#include <QtCore/QDir>


//OpenMS
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPASBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>

//STL
#include <iostream>
#include <map>
#include <vector>

#ifdef OPENMS_WINDOWSPLATFORM
#   ifndef _WIN32_WINNT
#       define _WIN32_WINNT 0x0501 // Win XP (and above)
#   endif
#   include <Windows.h>
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "TOPPAS";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage(Logger::LogStream& stream = OpenMS_Log_info)
{
  stream << "\n"
         << tool_name << " -- An assistant for GUI-driven TOPP workflow design." << "\n"
         << "\n"
         << "Usage:" << "\n"
         << " " << tool_name << " [options] [.toppas files]" << "\n"
         << "\n"
         << "Options are:" << "\n"
         << "  --help           Shows this help" << "\n"
         << "  --debug          Enables debug messages\n"
         << "  -ini <File>      Sets the INI file (default: ~/.TOPPAS.ini)" << "\n"
         << endl;
}

int main(int argc, const char** argv)
{
#ifdef OPENMS_WINDOWSPLATFORM
  qputenv("QT_QPA_PLATFORM", "windows:darkmode=0"); // disable dark mode on Windows, since our buttons etc are not designed for it
#endif

  // list of all the valid options
  std::map<std::string, std::string> valid_options, valid_flags, option_lists;
  valid_flags["--help"] = "help";
  valid_flags["--debug"] = "debug";
  valid_options["-ini"] = "ini";
  // invalid, but keep for now in order to inform users where to find this functionality now
  valid_options["-execute"] = "execute";
  valid_options["-out_dir"] = "out_dir";

  Param param;
  param.parseCommandLine(argc, argv, valid_options, valid_flags, option_lists);

  // '--help' given
  if (param.exists("help"))
  {
    print_usage();
    return 0;
  }

  // '-debug' given
  if (param.exists("debug"))
  {
    OPENMS_LOG_INFO << "Debug flag provided. Enabling 'OPENMS_LOG_DEBUG' ..." << std::endl;
    OpenMS_Log_debug.insert(cout); // allows to use OPENMS_LOG_DEBUG << "something" << std::endl;
  }

  // test if unknown options were given
  if (param.exists("unknown"))
  {
    // if TOPPAS is packed as Mac OS X bundle it will get a -psn_.. parameter by default from the OS
    // if this is the only unknown option it will be ignored .. maybe this should be solved directly
    // in Param.h
    if (!(String(param.getValue("unknown").toString()).hasSubstring("-psn") && !String(param.getValue("unknown").toString()).hasSubstring(", ")))
    {
      OPENMS_LOG_ERROR << "Unknown option(s) '" << param.getValue("unknown").toString() << "' given. Aborting!" << endl;
      print_usage(OpenMS_Log_error);
      return 1;
    }
  }

  try
  {

    if (param.exists("execute") || param.exists("out_dir"))
    {
      OPENMS_LOG_ERROR << "The parameters '-execute' and '-out_dir' are not valid anymore. This functionality has been moved to the ExecutePipeline tool." << endl;
      return 1;
    }

    QApplicationTOPP a(argc, const_cast<char**>(argv));
    a.connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));

    TOPPASBase mw;
    mw.show();

    a.connect(&a, &QApplicationTOPP::fileOpen, &mw, &TOPPASBase::openToppasFile);

    // Create the splashscreen that is displayed while the application loads (version is drawn dynamically)
    QPixmap qpm(":/TOPPAS_Splashscreen.png");
    QPainter pt_ver(&qpm);
    pt_ver.setFont(QFont("Helvetica [Cronyx]", 15, 2, true));
    pt_ver.setPen(Qt::black);
    // draw version number dynamcially on top left corner
    pt_ver.drawText(5, 5+15, VersionInfo::getVersion().toQString());
    QSplashScreen splash_screen(qpm);
    splash_screen.show();
    
    QApplication::processEvents();
    StopWatch stop_watch;
    stop_watch.start();

    if (param.exists("ini"))
    {
      mw.loadPreferences(param.getValue("ini").toString());
    }

    if (param.exists("misc"))
    {
      mw.loadFiles(ListUtils::toStringList<std::string>(param.getValue("misc")), &splash_screen);
    }
    else 
    {
      mw.newPipeline();
    }

    // We are about to show the application.
    // Proper time to remove the splashscreen, if at least 3 seconds have passed...
    while (stop_watch.getClockTime() < 3.0) /*wait*/

    {
    }
    stop_watch.stop();
    splash_screen.close();

#ifdef OPENMS_WINDOWSPLATFORM
    FreeConsole(); // get rid of console window at this point (we will not see any console output from this point on)
    AttachConsole(-1); // if the parent is a console, reattach to it - so we can see debug output - a normal user will usually not use cmd.exe to start a GUI)
#endif

    int result = a.exec();
    return result;
  }
  //######################## ERROR HANDLING #################################
  catch (Exception::UnableToCreateFile& e)
  {
    cout << String("Error: Unable to write file (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::FileNotFound& e)
  {
    cout << String("Error: File not found (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::FileNotReadable& e)
  {
    cout << String("Error: File not readable (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::FileEmpty& e)
  {
    cout << String("Error: File empty (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::ParseError& e)
  {
    cout << String("Error: Unable to read file (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::InvalidValue& e)
  {
    cout << String("Error: Invalid value (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::BaseException& e)
  {
    cout << String("Error: Unexpected error (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }

  return 1;
}
