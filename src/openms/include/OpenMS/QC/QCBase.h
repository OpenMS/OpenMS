// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <iostream>
#include <map>

namespace OpenMS
{
  class MSExperiment;

  /**
   * @brief This class serves as an abstract base class for all QC classes.
   *
   * It contains the important feature of encoding the input requirements
   * for a certain QC.
   */
  class OPENMS_DLLAPI QCBase
  {
  public:
    /**
     * @brief Enum to encode a file type as a bit.
     */
    enum class Requires 
      : UInt64 // 64 bit unsigned type for bitwise and/or operations (see below)
    {
      FAIL,         //< default, does not encode for anything
      RAWMZML,      //< mzML file is required
      POSTFDRFEAT,  //< Features with FDR-filtered pepIDs
      PREFDRFEAT,   //< Features with unfiltered pepIDs
      CONTAMINANTS, //< Contaminant Database
      TRAFOALIGN,   //< transformationXMLs for RT-alignment
      SIZE_OF_REQUIRES
    };
    /// strings corresponding to enum Requires
    static const std::string names_of_requires[];

    /**
     * @brief Map to find a spectrum via its NativeID
    */
    class OPENMS_DLLAPI SpectraMap
    {
    public:
      /// Constructor
      SpectraMap() = default;

      /// CTor which allows immediate indexing of an MSExperiment
      SpectraMap(const MSExperiment& exp);

      /// Destructor
      ~SpectraMap() = default;

      /// calculate a new map, delete the old one
      void calculateMap(const MSExperiment& exp);

      /// get index from identifier
      /// @throws Exception::ElementNotFound if @p identifier is unknown
      UInt64 at(const String& identifier) const;

      /// clear the map
      void clear();
      
      /// check if empty
      bool empty() const;
      
      /// get size of map
      Size size() const;

    private:
      std::map<String, UInt64> nativeid_to_index_; //< nativeID to index
    };


    /**
     @brief Storing a status of available/needed inputs (i.e. a set of Requires) as UInt64
    
     Conversion from a Requires enum is computed as `pow(2, r)`.
     Multiple Requires attributes can be computed by bitwise 'or'.

     Only allows assignment and bit operations with itself and an object
     of type Requires, i.e. not with any numeric types.

    **/
    class Status
    {
    public:
      /// stream output for Status
      friend std::ostream& operator<<(std::ostream& os, const Status& stat);

      /// Constructors
      Status() : value_(0)
      {}

      Status(const Requires& req)
      {
        value_ = getPow_(req);
      }

      Status(const Status& stat)
      {
        value_ = stat.value_;
      }

      /// Assignment
      Status& operator=(const Requires& req)
      {
        value_ = getPow_(req);
        return *this;
      }

      /// Destructor (default)
      ~Status() = default;

      // Equal
      bool operator==(const Status& stat) const
      {
        return (value_ == stat.value_);
      }

      Status& operator=(const Status& stat) = default;
      // Bitwise operators

      Status operator&(const Requires& req) const
      {
        Status s = *this;
        s.value_ &= getPow_(req);
        return s;
      }

      Status operator&(const Status& stat) const
      {
        Status s = *this;
        s.value_ &= stat.value_;
        return s;
      }

      Status& operator&=(const Requires& req)
      {
        value_ &= getPow_(req);
        return *this;
      }

      Status& operator&=(const Status& stat)
      {
        value_ &= stat.value_;
        return *this;
      }

      Status operator|(const Requires& req) const
      {
        Status s = *this;
        s.value_ |= getPow_(req);
        return s;
      }

      Status operator|(const Status& stat) const
      {
        Status s = *this;
        s.value_ |= stat.value_;
        return s;
      }

      Status& operator|=(const Requires& req)
      {
        value_ |= getPow_(req);
        return *this;
      }

      Status& operator|=(const Status& stat)
      {
        value_ |= stat.value_;
        return *this;
      }

      /**
       * @brief Check if input status fulfills requirement status.
       */
      bool isSuperSetOf(const Status& stat) const
      {
        return ((value_ & stat.value_) == stat.value_);
      }

    private:
      /// computes pow(2, r)
      UInt64 getPow_(const Requires& r) const
      {
        return UInt64(1) << UInt64 (r);
      }
      UInt64 value_;
    };

    /**
    * @brief Returns the name of the metric
    */
    virtual const String& getName() const = 0;
    
    /**
     *@brief Returns the input data requirements of the compute(...) function
     */
    virtual Status requires() const = 0;
    

    /**
     * @brief function, which iterates through all PeptideIdentifications of a given FeatureMap and applies a given lambda function
     *
     * The Lambda may or may not change the PeptideIdentification
     */

    template <typename MAP, typename T>
    static void iterateFeatureMap(MAP& fmap, T lambda)
    {
      for (auto& pep_id : fmap.getUnassignedPeptideIdentifications())
      {
        lambda(pep_id);
      }

      for (auto& features : fmap)
      {
        for (auto& pep_id : features.getPeptideIdentifications())
        {
          lambda(pep_id);
        }
      }
    }
  };

  inline std::ostream& operator<<(std::ostream& os, const QCBase::Status& stat)
  {
    return os << stat.value_;
  }
}
