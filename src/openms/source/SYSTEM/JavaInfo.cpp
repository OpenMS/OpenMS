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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/JavaInfo.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QProcess>
#include <QtCore/QDir>

namespace OpenMS
{

  bool JavaInfo::canRun(const String& java_executable, bool verbose_on_error)
  {
    QProcess qp;
    qp.start(java_executable.toQString(), QStringList() << "-version", QIODevice::ReadOnly);
    bool success = qp.waitForFinished();
    if (!success && verbose_on_error)
    {
        OPENMS_LOG_ERROR << "Java-Check:\n";
        if (qp.error() == QProcess::Timedout)
        {
          OPENMS_LOG_ERROR
            << "  Java was found at '" << java_executable << "' but the process timed out (can happen on very busy systems).\n"
            << "  Please free some resources or if you want to run the TOPP tool nevertheless set the TOPP tools 'force' flag in order to avoid this check." << std::endl;
        }
        else if (qp.error() == QProcess::FailedToStart)
        {
          OPENMS_LOG_ERROR
            << "  Java not found at '" << java_executable << "'!\n"
            << "  Make sure Java is installed and this location is correct.\n";
          if (QDir::isRelativePath(java_executable.toQString()))
          {
            static String path;
            if (path.empty())
            {
              path = getenv("PATH");
            }
            OPENMS_LOG_ERROR << "  You might need to add the Java binary to your PATH variable\n"
              << "  or use an absolute path+filename pointing to Java.\n" 
              << "  The current SYSTEM PATH is: '" << path << "'.\n\n"
  #ifdef __APPLE__
              << "  On MacOSX, application bundles change the system PATH; Open your executable (e.g. KNIME/TOPPAS/TOPPView) from within the bundle (e.g. ./TOPPAS.app/Contents/MacOS/TOPPAS) to preserve the system PATH or use an absolute path to Java!\n"
  #endif
              << std::endl;
          }
          else 
          {
            OPENMS_LOG_ERROR << "  You gave an absolute path to Java. Please check if it's correct.\n"
              << "  You can also try 'java' if your system path is correctly configured.\n"
              << std::endl;
          }
        }
        else
        {
          OPENMS_LOG_ERROR << "  Error executing '" << java_executable << "'!\n"
                    << "  Error description: '" << qp.errorString().toStdString() << "'.\n";
        }
    }
    return success;
  }

} // namespace OpenMS
