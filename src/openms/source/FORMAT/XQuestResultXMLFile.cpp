// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XQuestResultXMLHandler.h>
#include <OpenMS/ANALYSIS/XLMS/OpenProXLUtils.h>

namespace OpenMS
{
  XQuestResultXMLFile::XQuestResultXMLFile() :
    XMLFile("/SCHEMAS/xQuest_1_0.xsd", "1.0"),
    n_hits_(-1),
    cum_hits_(NULL)
  {
  }
  XQuestResultXMLFile::~XQuestResultXMLFile()
  {
    delete this->cum_hits_;
  }

  void XQuestResultXMLFile::load(const String & filename, std::vector< XQuestResultMeta > & metas,
                                 std::vector< std::vector < CrossLinkSpectrumMatch > > & csms,
                                 bool calc_cum_hits)
  {
   for(std::vector< XQuestResultMeta >::iterator it = metas.begin(); it != metas.end(); it++)
   {
     it->clearMetaInfo();
   }
   this->n_hits_ = 0;
   delete this->cum_hits_;

   if (calc_cum_hits)
   {
     this->cum_hits_ = new std::vector< int >;
   }
   else
   {
     this->cum_hits_ = NULL;
   }

   Internal::XQuestResultXMLHandler handler(filename, metas, csms, this->n_hits_, this->cum_hits_);
   this->parse_(filename, &handler);
  }

  void XQuestResultXMLFile::delete_cum_hits()
  {
    delete this->cum_hits_;
    this->cum_hits_ = NULL;
  }

  int XQuestResultXMLFile::get_n_hits() const
  {
    return this->n_hits_;
  }

  std::vector< int > * XQuestResultXMLFile::get_cum_hits() const
  {
    return this->cum_hits_;
  }
}
