// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer:$
// $Author: Adam Tenderholt $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/SYSTEM/PSCPSecureFileTransfer.h>

using namespace std;

namespace OpenMS
{
  PSCPSecureFileTransfer::PSCPSecureFileTransfer()
    : AbstractSecureFileTransfer()
  {
    process_.setProcessChannelMode(QProcess::ForwardedChannels);
  }

  PSCPSecureFileTransfer::PSCPSecureFileTransfer(String hostname, String username, String password)
    : AbstractSecureFileTransfer(hostname, username, password)
  {
    process_.setProcessChannelMode(QProcess::ForwardedChannels);
  }

  PSCPSecureFileTransfer::~PSCPSecureFileTransfer()
  {
  }

  bool PSCPSecureFileTransfer::downloadFile(String fromFilename, String toFilename)
  {
    QStringList arguments;
    arguments << "-l" << username_.toQString() << "-pw" << password_.toQString();
    arguments << hostname_.toQString() + ":" + fromFilename.toQString() << toFilename.toQString();

    process_.start("pscp", arguments);

    if (!process_.waitForStarted())
    {
      displayPSCPError_((int) process_.error());
      return false;
    }

    if (!process_.waitForFinished(-1))
    {
      displayPSCPError_((int) process_.error());
      return false;
    }

    return true;

  }

  bool PSCPSecureFileTransfer::uploadFile(String fromFilename, String toFileName)
  {

    QStringList arguments;
    arguments << "-l" << username_.toQString() << "-pw" << password_.toQString();
    arguments << fromFilename.toQString() << hostname_.toQString() + ":" + toFileName.toQString();

    process_.start("pscp", arguments);

    if (!process_.waitForStarted())
    {
      displayPSCPError_((int) process_.error());
      return false;
    }

    if (!process_.waitForFinished(-1))
    {
      displayPSCPError_((int) process_.error());
      return false;
    }

    return true;

  }

  void PSCPSecureFileTransfer::displayPSCPError_(int error)
  {
    LOG_ERROR << "\n***********************************************************************\n";
    LOG_ERROR << "Problem with PSCP process.\n\n";
    switch(error)
    {
    case QProcess::FailedToStart:
      LOG_ERROR << "The PSCP executable is missing from your path. Please download it from\n";
      LOG_ERROR << "http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html, and either\n";
      LOG_ERROR << "add its directory to your PATH or copy it to the same directory from which\n";
      LOG_ERROR << "you call the PeakInvestigator tool.";
      break;
    case QProcess::Crashed:
      LOG_ERROR << "The PSCP program crashed. Please consult Veritomyx for support.";
      break;
    case QProcess::Timedout:
    case QProcess::UnknownError:
      LOG_ERROR << "There is an unknown or timeout error for starting the PSCP program.\n";
      LOG_ERROR << "Please contact Veritomyx for support.";
      break;
    case QProcess::ReadError:
    case QProcess::WriteError:
      LOG_ERROR << "There was a problem reading/writing to the PSCP program. Please contact\n";
      LOG_ERROR << "Veritomyx for support.";
      break;
    }
    LOG_ERROR << "\n***********************************************************************\n" << endl;
  }

} //OpenMS
