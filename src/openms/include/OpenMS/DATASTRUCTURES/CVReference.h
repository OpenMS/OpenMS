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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_CVREFERENCE_H
#define OPENMS_DATASTRUCTURES_CVREFERENCE_H

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS
{
  /**
      @brief Controlled Vocabulary Reference

      Reference to a controlled vocabulary, defined in the first section of a mapping file.

      @ingroup Datastructures
  */

  class OPENMS_DLLAPI CVReference
  {
public:

    /// Default constructor
    CVReference();

    /// Copy constructor
    CVReference(const CVReference & rhs);

    /// Destructor
    virtual ~CVReference();

    /// Assignment operator
    CVReference & operator=(const CVReference & rhs);

    /** @name Accessors
    */
    //@{
    /// sets the name of the CV reference
    void setName(const String & name);

    /// returns the name of the CV reference
    const String & getName() const;

    /// sets the CV identifier which is referenced
    void setIdentifier(const String & identifier);

    /// returns the CV identifier which is referenced
    const String & getIdentifier() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const CVReference & rhs) const;

    /// inequality operator
    bool operator!=(const CVReference & rhs) const;
    //@}


protected:

    String name_;

    String identifier_;
  };


} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CVREFERENCE_H
