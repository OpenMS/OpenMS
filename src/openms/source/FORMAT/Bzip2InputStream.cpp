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


#include <OpenMS/FORMAT/Bzip2InputStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace xercesc;

namespace OpenMS
{
  Bzip2InputStream::Bzip2InputStream(const   String & file_name) :
    bzip2_(new Bzip2Ifstream(file_name.c_str())), file_current_index_(0)
  {
  }

  Bzip2InputStream::Bzip2InputStream(const   char * file_name) :
    bzip2_(new Bzip2Ifstream(file_name)), file_current_index_(0)
  {
  }

/*	Bzip2InputStream::Bzip2InputStream()
    :bzip2_(NULL)
    {

    }*/

  Bzip2InputStream::~Bzip2InputStream()
  {
    delete bzip2_;
  }

  XMLSize_t Bzip2InputStream::readBytes(XMLByte * const  to_fill, const XMLSize_t  max_to_read)
  {
    //  Figure out whether we can really read.
    if (bzip2_->streamEnd())
    {
      return 0;
    }

    unsigned char * fill_it = static_cast<unsigned char *>(to_fill);
    XMLSize_t actual_read = (XMLSize_t) bzip2_->read((char *)fill_it, static_cast<const size_t>(max_to_read));
    file_current_index_ += actual_read;
    return actual_read;
  }

  const XMLCh * Bzip2InputStream::getContentType() const
  {
    return nullptr;
  }

} // namespace OpenMS
