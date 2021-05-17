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

#pragma once

#include <OpenMS/APPLICATIONS/TOPPBase.h>


namespace OpenMS
{
  
  /**
    @brief Base class for Search Engine Adapters

    It is build on top of TOPPBase and provides convenience functions for regular tasks in SearchEngines.

    This base class enforces a common parameter scheme upon each adapter.
    E.g. '-database' and '-in'.
    This might be extended/changed in the future.
    
  */
  class OPENMS_DLLAPI SearchEngineBase : public TOPPBase
  {
  public:
    /// No default constructor
    SearchEngineBase() = delete;

    /// No default copy constructor.
    SearchEngineBase(const SearchEngineBase&) = delete;

    /**
      @brief Constructor

      Must match TOPPBase' Ctor!

      @param name Tool name.
      @param description Short description of the tool (one line).
      @param official If this is an official TOPP tool contained in the OpenMS/TOPP release.
      If @em true the tool name is checked against the list of TOPP tools and a warning printed if missing.

      @param citations Add one or more citations if they are associated specifically to this TOPP tool; they will be printed during --help
    */
    SearchEngineBase(const String& name, const String& description, bool official = true, const std::vector<Citation>& citations = {}, bool toolhandler_test = true);

    /// Destructor
    virtual ~SearchEngineBase();

    /**
      @brief Reads the '-in' argument from internal parameters (usually an mzML file) and checks if MS2 spectra are present and are centroided.

      If the file is an mzML file, the spectra annotation can be checked. If no MS2 or profile MS2 data is found, an exception is thrown.
      If the file is any other format, the overhead of reading in the file is too large and we just issue a general warning that centroided data should be used.

      @param ms_level The MS level to check for their type (centroided/profile)

      @return A filename (might be a relative or absolute path)

      @throws OpenMS::Exception::FileEmpty if no spectra are found (mzML only)
      @throws OpenMS::Exception::IllegalArgument if spectra are not centroided (mzML only)

    */
    String getRawfileName(int ms_level = 2) const;


    /**
      @brief Reads the '-database' argument from internal parameters (or from @p db) and tries to find the db in search directories (if it cannot be found immediately). If not found, an exception is thrown.
      
      @param db [Optional] Instead of reading the '-database', you can provide a custom name here (might be required for special db formats, see OMSSA)
      @return filename for DB (might be a relative or absolute path)
 
      @throws OpenMS::Exception::FileNotFound if database name could not be resolved
    */
    String getDBFilename(String db = "") const;
  }; // end SearchEngineBase

}   // end NS OpenMS
