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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_TRANSFORMATIONXMLFILE_H
#define OPENMS_FORMAT_TRANSFORMATIONXMLFILE_H

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Used to load and store TransformationXML files

    This class is used to load and store documents that implement the schema of
    TransformationXML files.

    A documented schema for this format can be found at http://open-ms.sourceforge.net/schemas/.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI TransformationXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// Constructor
    TransformationXMLFile();

    /**
    @brief Loads the transformation from an TransformationXML file

    The information is read in and the information is stored in the
    corresponding variables

    @exception Exception::FileNotFound is thrown if the file could not be opened
    @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, TransformationDescription& transformation, bool fit_model=true);

    /**
    @brief Stores the data in an TransformationXML file

    The data is read in and stored in the file named 'filename'.

    @exception Exception::UnableToCreateFile is thrown if the file could not be created
    @exception Exception::IllegalArgument is thrown if unsupported parameter types have been set
    */
    void store(String filename, const TransformationDescription& transformation);

protected:
    // Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

    /// @name Members for use during loading data
    //@{
    /// Param to fill in during parse
    Param params_;
    /// Data vector
    TransformationDescription::DataPoints data_;
    /// Model type
    String model_type_;
    //@}

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_TRANSFORMATIONXMLFILE_H
