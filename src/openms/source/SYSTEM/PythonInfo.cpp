// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/PythonInfo.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QProcess>
#include <QtCore/QDir>

#include <sstream>

using namespace std;

namespace OpenMS
{
  bool PythonInfo::canRun(String& python_executable, String& error_msg)
  {
    stringstream ss;
    String py_original = python_executable;
    if (!File::findExecutable(python_executable))
    {
      ss << "  Python not found at '" << python_executable << "'!\n"
         << "  Make sure Python is installed and this location is correct.\n";
      if (QDir::isRelativePath(python_executable.toQString()))
      {
        static String path;
        if (path.empty())
        {
          path = getenv("PATH");
        }
        ss << "  You might need to add the Python binary to your PATH variable\n"
           << "  or use an absolute path+filename pointing to Python.\n"
           << "  The current SYSTEM PATH is: '" << path << "'.\n\n";
#ifdef __APPLE__
        ss << "  On MacOSX, application bundles change the system PATH; Open your executable (e.g. KNIME/TOPPAS/TOPPView) from within the bundle (e.g. ./TOPPAS.app/Contents/MacOS/TOPPAS) to preserve the system PATH or use an absolute path to Python!\n";
#endif
      }
      error_msg = ss.str();
      return false;
    }

    if (python_executable != py_original)
    {
      ss << "Python executable ('" << py_original << "') resolved to '" << python_executable << "'\n";
    }

    QProcess qp;
    qp.start(python_executable.toQString(), QStringList() << "--version", QIODevice::ReadOnly);
    bool success = qp.waitForFinished();
    if (!success)
    {
      if (qp.error() == QProcess::Timedout)
      {
        ss << "  Python was found at '" << python_executable << "' but the process timed out (can happen on very busy systems).\n"
           << "  Please free some resources or if you want to run the TOPP tool nevertheless set the TOPP tools 'force' flag in order to avoid this check.\n";
      }
      else if (qp.error() == QProcess::FailedToStart)
      {
        ss << "  Python found at '" << python_executable << "' but failed to run!\n"
           << "  Make sure you have the rights to execute this binary file.\n";
      }
      else
      {
        ss << "  Error executing '" << python_executable << "'!\n"
           << "  Error description: '" << qp.errorString().toStdString() << "'.\n";
      }
    }

    error_msg = ss.str();
    return success;
  }

  bool PythonInfo::isPackageInstalled(const String& python_executable, const String& package_name)
  {
    QProcess qp;
    qp.start(python_executable.toQString(), QStringList() << "-c" << (String("import ") + package_name).c_str(), QIODevice::ReadOnly);
    bool success = qp.waitForFinished();
    return (success && qp.exitStatus() == QProcess::ExitStatus::NormalExit && qp.exitCode() == 0);
  }

  String PythonInfo::getVersion(const String& python_executable)
  {
    String v;
    QProcess qp;
    qp.start(python_executable.toQString(), QStringList() << "--version", QIODevice::ReadOnly);
    bool success = qp.waitForFinished();
    if (success && qp.exitStatus() == QProcess::ExitStatus::NormalExit && qp.exitCode() == 0)
    {
      v = qp.readAllStandardOutput().toStdString(); // some pythons report is on stdout
      v += qp.readAllStandardError().toStdString();  // ... some on stderr
      v.trim(); // remove '\n'
    }
    return v;
  }

} // namespace OpenMS
