// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/CompressedInputSource.h>

#include <OpenMS/FORMAT/GzipInputStream.h>
#include <OpenMS/FORMAT/Bzip2InputStream.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <xercesc/internal/MemoryManagerImpl.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLUniDefs.hpp>


using namespace xercesc;
namespace OpenMS
{

  CompressedInputSource::CompressedInputSource(const String & file_path, const String & header, MemoryManager * const manager) :
    xercesc::InputSource(manager),
    head_(header)
  {
    if (head_.size() < 2)
    {
      head_ = "\0\0";
    }
    //
    //  If the path is relative, then complete it according to the current
    //  working directory rules of the current platform. Else, just take
    //  it as is.
    //
    Internal::StringManager strman;
    auto file = strman.convert(file_path.c_str());
    if (xercesc::XMLPlatformUtils::isRelative(file.c_str(), manager))
    {
      XMLCh * curDir = xercesc::XMLPlatformUtils::getCurrentDirectory(manager);

      XMLSize_t curDirLen = XMLString::stringLen(curDir);
      XMLSize_t filePathLen = XMLString::stringLen(file.c_str());
      XMLCh * fullDir = (XMLCh *) manager->allocate
                        (
        (curDirLen + filePathLen + 2) * sizeof(XMLCh)
                        ); //new XMLCh [ curDirLen + filePathLen + 2];

      XMLString::copyString(fullDir, curDir);
      fullDir[curDirLen] = chForwardSlash;
      XMLString::copyString(&fullDir[curDirLen + 1], file.c_str());

      XMLPlatformUtils::removeDotSlash(fullDir, manager);
      XMLPlatformUtils::removeDotDotSlash(fullDir, manager);

      setSystemId(fullDir);

      manager->deallocate(curDir);  //delete [] curDir;
      manager->deallocate(fullDir);  //delete [] fullDir;
    }
    else
    {
      XMLCh * tmpBuf = XMLString::replicate(file.c_str(), manager);
      XMLPlatformUtils::removeDotSlash(tmpBuf, manager);
      setSystemId(tmpBuf);
      manager->deallocate(tmpBuf);  //delete [] tmpBuf;
    }
  }

  CompressedInputSource::CompressedInputSource(const XMLCh * const file, const String & header, MemoryManager * const manager) :
    xercesc::InputSource(manager),
    head_(header)
  {
    if (head_.size() < 2)
    {
      head_ = "\0\0";
    }
    //
    //  If the path is relative, then complete it according to the current
    //  working directory rules of the current platform. Else, just take
    //  it as is.
    //
    if (xercesc::XMLPlatformUtils::isRelative(file, manager))
    {
      XMLCh * curDir = xercesc::XMLPlatformUtils::getCurrentDirectory(manager);

      XMLSize_t curDirLen = XMLString::stringLen(curDir);
      XMLSize_t filePathLen = XMLString::stringLen(file);
      XMLCh * fullDir = (XMLCh *) manager->allocate
                        (
        (curDirLen + filePathLen + 2) * sizeof(XMLCh)
                        ); //new XMLCh [ curDirLen + filePathLen + 2];

      XMLString::copyString(fullDir, curDir);
      fullDir[curDirLen] = chForwardSlash;
      XMLString::copyString(&fullDir[curDirLen + 1], file);

      XMLPlatformUtils::removeDotSlash(fullDir, manager);
      XMLPlatformUtils::removeDotDotSlash(fullDir, manager);

      setSystemId(fullDir);

      manager->deallocate(curDir);  //delete [] curDir;
      manager->deallocate(fullDir);  //delete [] fullDir;
    }
    else
    {
      XMLCh * tmpBuf = XMLString::replicate(file, manager);
      XMLPlatformUtils::removeDotSlash(tmpBuf, manager);
      setSystemId(tmpBuf);
      manager->deallocate(tmpBuf);  //delete [] tmpBuf;
    }
  }

  CompressedInputSource::~CompressedInputSource()
  {
  }

  BinInputStream * CompressedInputSource::makeStream() const
  {
    if (head_[0] == 'B' && head_[1] == 'Z')
    {
      Bzip2InputStream * retStrm = new Bzip2InputStream(Internal::StringManager().convert(getSystemId()));
      if (!retStrm->getIsOpen())
      {
        delete retStrm;
        return 0;
      }
      return retStrm;
    }
    else /*     (bz[0] == g1 && bz[1] == g2), where char g1 = 0x1f and char g2 = 0x8b */
    {
      GzipInputStream * retStrm = new GzipInputStream(Internal::StringManager().convert(getSystemId()));
      if (!retStrm->getIsOpen())
      {
        delete retStrm;
        return 0;
      }
      return retStrm;
    }

  }

} // namespace OpenMS
