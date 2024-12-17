// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

/**
  @page TOPP_TOPPView TOPPView

  TOPPView is a viewer for MS and HPLC-MS data. It can be used to inspect files in mzML, mzData, mzXML
  and several other file formats. It also supports viewing data from an %OpenMS database.
  The following figure shows two instances of TOPPView displaying a HPLC-MS map and a MS raw spectrum:

  @image html TOPPView.png

  More information about TOPPView can be found on the OpenMS ReadTheDocs
  page: https://openms.readthedocs.io/en/latest/openms-applications-and-tools/visualize-with-openms.html

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_TOPPView.cli
  
  Note: By default, TOPPView scans for novel TOPP tools if there has been a version update. To force a rescan you
  can pass the --force flag. To skip the scan for tools, you can pass the --skip_tool_scan flag.
*/

//QT
#include <QtWidgets/QSplashScreen>
#include <QMessageBox>

//OpenMS
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/SYSTEM/StopWatch.h>

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
const char* tool_name = "TOPPView";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage()
{
  cerr << endl
       << tool_name << " -- A viewer for mass spectrometry data." << "\n"
       << "\n"
       << "Usage:" << "\n"
       << " " << tool_name << " [options] [files]" << "\n"
       << "\n"
       << "Options are:" << "\n"
       << "  --help           Shows this help" << "\n"
       << "  -ini <File>      Sets the INI file (default: ~/.TOPPView.ini)" << "\n"
       << "  --force          Forces scan for new tools" << "\n"
       << "  --skip_tool_scan Skips scan for new tools" << "\n"
       << "\n"
       << "Hints:" << "\n"
       << " - To open several files in one window put a '+' in between the files." << "\n"
       << " - '@bw' after a map file displays the dots in a white to black gradient." << "\n"
       << " - '@bg' after a map file displays the dots in a grey to black gradient." << "\n"
       << " - '@b'  after a map file displays the dots in black." << "\n"
       << " - '@r'  after a map file displays the dots in red." << "\n"
       << " - '@g'  after a map file displays the dots in green." << "\n"
       << " - '@m'  after a map file displays the dots in magenta." << "\n"
       << " - Example: '" << tool_name << " 1.mzML + 2.mzML @bw + 3.mzML @bg'" << "\n"
       << endl;
}

int main(int argc, const char** argv)
{
 #ifdef OPENMS_WINDOWSPLATFORM
  qputenv("QT_QPA_PLATFORM", "windows:darkmode=0"); // disable dark mode on Windows, since our buttons etc are not designed for it
#endif

  //list of all the valid options
  std::map<std::string, std::string> valid_options, valid_flags, option_lists;
  valid_flags["--help"] = "help";
  valid_flags["--force"] = "force";
  valid_flags["--skip_tool_scan"] = "skip_tool_scan";
  valid_flags["--debug"] = "debug";
  valid_options["-ini"] = "ini";

  Param param;
  param.parseCommandLine(argc, argv, valid_options, valid_flags, option_lists);

  // '--help' given
  if (param.exists("help"))
  {
    print_usage();
    return 0;
  }

  // test if unknown options were given
  if (param.exists("unknown"))
  {
    // if TOPPView is packed as Mac OS X bundle it will get a -psn_.. parameter by default from the OS
    // if this is the only unknown option it will be ignored .. maybe this should be solved directly
    // in Param.h
    if (!(String(param.getValue("unknown").toString()).hasSubstring("-psn") && !String(param.getValue("unknown").toString()).hasSubstring(", ")))
    {
      cout << "Unknown option(s) '" << param.getValue("unknown").toString() << "' given. Aborting!" << endl;
      print_usage();
      return 1;
    }
  }

  try
  {

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

    QApplicationTOPP a(argc, const_cast<char**>(argv));
    a.connect(&a, &QApplicationTOPP::lastWindowClosed, &a, &QApplicationTOPP::quit);

    TOPPViewBase::TOOL_SCAN mode = TOPPViewBase::TOOL_SCAN::SCAN_IF_NEWER_VERSION;
    if (param.exists("force"))
    {
      mode = TOPPViewBase::TOOL_SCAN::FORCE_SCAN;
    }
    else if (param.exists("skip_tool_scan"))
    {
      mode = TOPPViewBase::TOOL_SCAN::SKIP_SCAN;
    }

    TOPPViewBase::VERBOSITY verbosity = TOPPViewBase::VERBOSITY::DEFAULT;
    if (param.exists("debug"))
    {
      verbosity = TOPPViewBase::VERBOSITY::VERBOSE;
    }

    TOPPViewBase tb(mode, verbosity);
    a.connect(&a, &QApplicationTOPP::fileOpen, &tb, &TOPPViewBase::openFile);
    tb.show();

    // Create the splashscreen that is displayed while the application loads (version is drawn dynamically)
    QPixmap qpm(":/TOPPView_Splashscreen.png");
    QPainter pt_ver(&qpm);
    pt_ver.setFont(QFont("Helvetica [Cronyx]", 15, 2, true));
    pt_ver.setPen(Qt::black);
    // draw version number dynamcially on top left corner
    pt_ver.drawText(5, 5 + 15, VersionInfo::getVersion().toQString());
    QSplashScreen splash_screen(qpm);
    splash_screen.show();
    
    QApplication::processEvents();
    StopWatch stop_watch;
    stop_watch.start();

    if (param.exists("ini"))
    {
      tb.loadPreferences(param.getValue("ini").toString());
    }

    //load command line files
    if (param.exists("misc"))
    {
      tb.loadFiles(ListUtils::toStringList<std::string>(param.getValue("misc")), &splash_screen);
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
    return a.exec();
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
