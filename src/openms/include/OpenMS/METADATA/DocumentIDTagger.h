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
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Tags OpenMS file containers with a DocumentID

    Intended usage is from within a TOPP tool. An instance of this class
        is present in TOPPBase and can be used by all derived TOPP tools
        to assign a unique ID which is fetched from an ID pool in ./share/OpenMS/IDPool/.

        @ingroup Metadata
  */
  class OPENMS_DLLAPI DocumentIDTagger
  {
private:
    /// Constructor (declared away)
    DocumentIDTagger();

    /**
        @brief retrieve an ID from the pool

        Uses boost file locks to safely retrieve an ID from an ID pool.

        @param id Unique identifier returned from ID pool
        @param free Number of available identifiers in ID pool (before this query)
        @param idcount_only Only count available identifiers, do NOT retrieve one
                     (the id string will nevertheless be filled)

        Return true if all file operations could be executed successfully (this does not imply there was an ID left over - check free>0)
    */
    bool getID_(String & id, Int & free, bool idcount_only) const;

public:
    /// Constructor
    DocumentIDTagger(String toolname);
    /// Copy constructor
    DocumentIDTagger(const DocumentIDTagger & source) = default;
    /// Destructor
    ~DocumentIDTagger();

    /// Assignment operator
    DocumentIDTagger & operator=(const DocumentIDTagger & source) = default;

    /// Equality operator
    bool operator==(const DocumentIDTagger & source) const;
    /// Equality operator
    bool operator!=(const DocumentIDTagger & source) const;


    /**
        @brief Return the file used as ID pool

        The default ID pool file is in /share/OpenMS/IDPool/IDPool.txt
        A custom file can be set by setIDPoolFile()
    */
    String getPoolFile() const;

    /// Set the file used as ID pool
    void setPoolFile(const String & file);

    /**
        @brief Tags any structure which is derived from DocumentIdentifier with a unique tag

        Tags any structure which is derived from DocumentIdentifier with a unique tag
        Returns true if ID could be assigned, otherwise an Exception::DepletedIDPool is thrown

        @param map Some class (derived from a DocumentIdentifier class) which needs a unique id
        @exception Exception::DepletedIDPool when no identifier (for whatever reason) could be acquired
    */
    bool tag(DocumentIdentifier & map) const;

    /**
        @brief return the number of available IDs in the pool.

        Retrieve the number of available IDs in the pool.
        Returns true of count was successful, false otherwise (locking error, file creation error ...)

        @param free Number of available identifiers. You should worry if it's 0!
    */
    bool countFreeIDs(Int & free) const;

protected:
    /// name of the calling TOPP tool
    String toolname_;

    /// location of the ID pool
    String pool_file_;
  };

} // namespace OpenMS

