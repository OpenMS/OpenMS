// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Nico Pfeifer $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MASCOTXMLFILE_H
#define OPENMS_FORMAT_MASCOTXMLFILE_H

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/MascotXMLHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>


namespace OpenMS
{
  class ProteinIdentification;

  /**
    @brief Used to load MascotXML files

    This class is used to load documents that implement
    the schema of MascotXML files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI MascotXMLFile :
    public Internal::XMLFile
  {
public:

    typedef Internal::MascotXMLHandler::RTMapping RTMapping;

    /// Constructor
    MascotXMLFile();

    /**
      @brief loads data from a Mascot XML file

      @param filename the file to be loaded
      @param protein_identification protein identifications belonging to the whole experiment
      @param id_data the identifications with m/z and RT
      @param rt_mapping An optional mapping of scan indices to RT, in case the file only contains scan numbers
      @param scan_regex An optional regular expression used to extract the scan numbers

      @exception Exception::FileNotFound is thrown if the file does not exists.
      @exception Exception::ParseError is thrown if the file does not suit to the standard.

      This method serves to read in a Mascot XML file. The information can be
      retrieved via the load function.
    */
    void load(const String& filename,
              ProteinIdentification& protein_identification,
              std::vector<PeptideIdentification>& id_data,
              const RTMapping& rt_mapping = RTMapping(),
              const String& scan_regex = "");

    /**
      @brief loads data from a Mascot XML file

      @param filename the file to be loaded
      @param protein_identification protein identifications belonging to the whole experiment
      @param id_data the identifications with m/z and RT
      @param peptides a map of modified peptides identified by the String title
      @param rt_mapping An optional mapping of scan indices to RT, in case the file only contains scan numbers
      @param scan_regex An optional regular expression used to extract the scan numbers

      @exception Exception::FileNotFound is thrown if the file does not exists.
      @exception Exception::ParseError is thrown if the file does not suit to the standard.

      This method serves to read in a Mascot XML file. The information can be
      retrieved via the load function.
    */
    void load(const String& filename,
              ProteinIdentification& protein_identification,
              std::vector<PeptideIdentification>& id_data, 
              std::map<String, std::vector<AASequence> >& peptides, 
              const RTMapping& rt_mapping = RTMapping(),
              const String& scan_regex = "");

    /**
      @brief Generates a mapping between scan numbers and retention times in raw data

      @param begin Iterator to the first spectrum
      @param end Iterator past the last spectrum
      @param rt_mapping Output mapping

      The mapping can be used to infer retention times of identifications when a Mascot XML file is loaded.
    */  
    static void generateRTMapping(const MSExperiment<>::ConstIterator begin, 
                                  const MSExperiment<>::ConstIterator end, 
                                  RTMapping& rt_mapping);

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MASCOTXMLFILE_H
