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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>

namespace OpenMS
{
  namespace Internal
  {

    String MzMLHandlerHelper::getCompressionTerm_(const PeakFileOptions& opt, MSNumpressCoder::NumpressConfig np, bool use_numpress)
    {
      if (np.np_compression != MSNumpressCoder::NONE && opt.getCompression() )
      {
        // TODO check if zlib AND numpress are allowed at the same time by the standard ... 
        // It is technically possible but
        //
        // MUST supply a *child* term of MS:1000572 (binary data compression type) only once
        //
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Cannot have numpress and zlib compression at the same time", "numpress, zlib");
      }

      if (np.np_compression == MSNumpressCoder::NONE || ! use_numpress)
      {
        if (opt.getCompression())
        {
          return "<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\" />";
        }
        else
        {
          return "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" />";
        }
      }
      else if (np.np_compression == MSNumpressCoder::LINEAR)
      {
        return "<cvParam cvRef=\"MS\" accession=\"MS:1002312\" name=\"MS-Numpress linear prediction compression\" />";
      }
      else if (np.np_compression == MSNumpressCoder::PIC)
      {
        return "<cvParam cvRef=\"MS\" accession=\"MS:1002313\" name=\"MS-Numpress linear prediction compression\" />";
      }
      else if (np.np_compression == MSNumpressCoder::SLOF)
      {
        return "<cvParam cvRef=\"MS\" accession=\"MS:1002314\" name=\"MS-Numpress short logged float compression\" />";
      }
      else
      {
        // default
        return "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" />";
      }
    }

    void MzMLHandlerHelper::writeFooter_(std::ostream& os, const PeakFileOptions& options_, 
      std::vector< std::pair<std::string, long> > & spectra_offsets,
      std::vector< std::pair<std::string, long> > & chromatograms_offsets)
    {
      os << "\t</run>\n";
      os << "</mzML>";

      if (options_.getWriteIndex())
      {
        int indexlists = (int) !spectra_offsets.empty() + (int) !chromatograms_offsets.empty();

        long indexlistoffset = os.tellp();
        os << "\n";
        // NOTE: indexList is required, so we need to write one 
        os << "  <indexList count=\"" << indexlists << "\">\n";
        if (!spectra_offsets.empty())
        {
          os << "    <index name=\"spectrum\">\n";
          for (Size i = 0; i < spectra_offsets.size(); i++)
          {
            os << "      <offset idRef=\"" << spectra_offsets[i].first << "\">" << spectra_offsets[i].second << "</offset>\n";
          }
          os << "    </index>\n";
        }
        if (!chromatograms_offsets.empty())
        {
          os << "    <index name=\"chromatogram\">\n";
          for (Size i = 0; i < chromatograms_offsets.size(); i++)
          {
            os << "      <offset idRef=\"" << chromatograms_offsets[i].first << "\">" << chromatograms_offsets[i].second << "</offset>\n";
          }
          os << "    </index>\n";
        }
        if (indexlists == 0)
        {
          // dummy: at least one index subelement is required by the standard,
          // and at least one offset element is required so we need to handle
          // the case where no spectra/chromatograms are present.
          os << "    <index name=\"dummy\">\n";
            os << "      <offset idRef=\"dummy\">-1</offset>\n";
          os << "    </index>\n";
        }
        os << "  </indexList>\n";
        os << "  <indexListOffset>" << indexlistoffset << "</indexListOffset>\n";
        os << "<fileChecksum>";

        // TODO calculate checksum here:
        //  SHA-1 checksum from beginning of file to end of 'fileChecksum' open tag.
        String sha1_checksum = "0";
        os << sha1_checksum << "</fileChecksum>\n";

        os << "</indexedmzML>";
      }
    }

  }
} // namespace OpenMS
