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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_INDEXEDMZMLDECODER_H
#define OPENMS_FORMAT_HANDLERS_INDEXEDMZMLDECODER_H

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS 
{

  /**
    @brief A class to analyze indexedmzML files and extract the offsets of individual tags

    Specifically, this class allows one to extract the offsets of the <indexList>
    tag and of all <spectrum> and <chromatogram> tag using the indices found at
    the end of the indexedmzML XML structure.

    While findIndexListOffset tries extracts the offset of the indexList tag from
    the last 1024 bytes of the file, this offset allows the function parseOffsets
    to extract all elements contained in the <indexList> tag and thus get access
    to all spectra and chromatogram offsets.

  */
  class OPENMS_DLLAPI IndexedMzMLDecoder
  {
  public:

    /// The vector containing binary offsets
    typedef std::vector< std::pair<std::string, long> > OffsetVector;

    /**
      @brief Tries to extract the offsets of all spectra and chromatograms from an indexedmzML.

      Given the start of the <indexList> element, this function tries to read
      this tag from the given the indexedmzML file. It stores the result in the
      spectra and chromatogram offset vectors.

      @param in Filename of the input indexedmzML file
      @param indexoffset Offset at which position in the file the XML tag "<indexList" is expected to occur
      @param spectra_offsets Output vector containing the positions of all spectra in the file
      @param chromatograms_offsets Output vector containing the positions of all chromatograms in the file

      @return 0 in case of success and -1 otherwise (failure, no offset was found)

    */
    int parseOffsets(String filename, int indexoffset, OffsetVector & spectra_offsets, OffsetVector& chromatograms_offsets);

    /**
      @brief Tries to extract the indexList offset from an indexedmzML.

      This function reads by default the last few (1024) bytes of the given
      input file and tries to read the content of the <indexListOffset> tag.
      The idea is that somewhere in the last parts of the file specified by the
      input string, the string <indexListOffset>xxx</indexListOffset> occurs.
      This function returns the xxx part converted to an integer.

      @note Since this function cannot determine where it will start reading
      the XML, no regular XML parser can be used for this. Therefore it uses
      regex to do its job. It matches the <indexListOffset> part and any
      numerical characters that follow. 

      @param in Filename of the input indexedmzML file
      @param buffersize How many bytes of the input file should be searched for the tag

      @return A positive integer containing the content of the indexListOffset tag, returns -1 in case of failure no tag was found

    */
    int findIndexListOffset(String filename, int buffersize = 1023);

  protected:

    /**
      @brief Extract data from a string containing an <indexList> tag

      This function parses the contained <offset> tags inside the indexList tag
      and stores the contents in the spectra and chromatogram offset vectors.

      This function expects an input string that contains a root XML tag and as
      one of its child an <indexList> tag as defined by the mzML 1.1.0 index
      wrapper schema. Usually the root would be an indexedmzML tag and _must_
      contain an indexList tag, while the dx:mzML, indexListOffset and
      fileChecksum are optional(their presence is not checked).

      Still this means, don't stick non-valid XML in here (e.g. non matching
      open/close tags).  Usually this means that you will at least have to add
      an opening </indexedmzML>. Valid input for this function would for
      example be:

      @code

      <indexedmzML>
        <indexList count="1">
          <index name="chromatogram">
            <offset idRef="1">9752</offset>
          </index>
        </indexList>
        <indexListOffset>26795</indexListOffset>
      <fileChecksum>0</fileChecksum>
      </indexedmzML>

      @endcode
      
      @param in String containing the XML with a indexedmzML parent and an indexList child tag
      @param spectra_offsets Output vector containing the positions of all spectra in the file
      @param chromatograms_offsets Output vector containing the positions of all chromatograms in the file
    */
    int domParseIndexedEnd(std::string in, OffsetVector & spectra_offsets, OffsetVector& chromatograms_offsets);
  };

}
#endif
