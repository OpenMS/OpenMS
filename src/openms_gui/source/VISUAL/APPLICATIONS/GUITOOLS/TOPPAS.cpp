// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
#include <QtCore/QDir>


//OpenMS
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/VISUAL/MISC/InteractiveSplashScreen.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPASBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>

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

void print_usage(Logger::LogStream& stream = getGlobalLogInfo())
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
         << "\n"
         << "Note: Qt command line options (e.g. '-style <style>' or '-stylesheet <file>')\n"
         << "      are supported as well and are passed on to Qt.\n"
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
  // (handled before constructing a QApplication, so that '--help' also works in headless environments)
  if (param.exists("help"))
  {
    print_usage();
    return 0;
  }

  try
  {

    QApplicationTOPP a(argc, const_cast<char**>(argv));
    a.connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));

    // Qt has now consumed (and removed from argc/argv) the command line arguments it recognizes,
    // e.g. '-style', '-stylesheet', '-platform', ... (see https://doc.qt.io/qt-5/qapplication.html#QApplication).
    // This allows users to customize the GUI appearance. We therefore re-parse the now reduced command
    // line and only afterwards check for unknown options, so that Qt arguments are not mistaken for them.
    param.clear();
    param.parseCommandLine(argc, argv, valid_options, valid_flags, option_lists);

    // '-debug' given
    if (param.exists("debug"))
    {
      OPENMS_LOG_INFO << "Debug flag provided. Enabling 'OPENMS_LOG_DEBUG' ..." << std::endl;
      getGlobalLogDebug().insert(cout); // allows to use OPENMS_LOG_DEBUG << "something" << std::endl;
    }

    // test if unknown options were given
    if (param.exists("unknown"))
    {
      // if TOPPAS is packed as Mac OS X bundle it will get a -psn_.. parameter by default from the OS
      // if this is the only unknown option it will be ignored .. maybe this should be solved directly
      // in Param.h
      if (!(StringUtils::hasSubstring(std::string(param.getValue("unknown").toString()), "-psn") && !StringUtils::hasSubstring(std::string(param.getValue("unknown").toString()), ", ")))
      {
        OPENMS_LOG_ERROR << "Unknown option(s) '" << param.getValue("unknown").toString() << "' given. Aborting!" << endl;
        print_usage(getGlobalLogError());
        return 1;
      }
    }

    if (param.exists("execute") || param.exists("out_dir"))
    {
      OPENMS_LOG_ERROR << "The parameters '-execute' and '-out_dir' are not valid anymore. This functionality has been moved to the ExecutePipeline tool." << endl;
      return 1;
    }

    TOPPASBase mw;
    mw.show();

    a.connect(&a, &QApplicationTOPP::fileOpen, &mw, &TOPPASBase::openToppasFile);

    // Create the splashscreen that is displayed while the application loads (version is drawn dynamically)
    QPixmap qpm(":/TOPPAS_Splashscreen.png");
    QPainter pt_ver(&qpm);
    pt_ver.setFont(QFont("Helvetica [Cronyx]", 15, 2, true));
    pt_ver.setPen(Qt::black);
    // draw version number dynamcially on top left corner
    pt_ver.drawText(5, 5+15, toQString(VersionInfo::getVersion()));
    InteractiveSplashScreen splash_screen(qpm);
    splash_screen.show();

    QApplication::processEvents();

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

    // Keep the splashscreen up for at least 3 seconds so it can be read, but let the user
    // dismiss it earlier with a mouse click or key press. The event loop stays responsive.
    splash_screen.showFor(3.0);


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
    cout << std::string("Error: Unable to write file (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::FileNotFound& e)
  {
    cout << std::string("Error: File not found (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::FileNotReadable& e)
  {
    cout << std::string("Error: File not readable (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::FileEmpty& e)
  {
    cout << std::string("Error: File empty (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::ParseError& e)
  {
    cout << std::string("Error: Unable to read file (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::InvalidValue& e)
  {
    cout << std::string("Error: Invalid value (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::BaseException& e)
  {
    cout << std::string("Error: Unexpected error (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }

  return 1;
}
