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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MODELDESCRIPTION_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MODELDESCRIPTION_H

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/CONCEPT/Factory.h>

#include <sstream>

namespace OpenMS
{
  /**
      @brief Stores the name and parameters of a model.

      This class also allows reconstruction of the model.

      @see BaseModel
  */
  template <UInt D>
  class ModelDescription
  {
public:

    /// Default constructor
    ModelDescription() :
      name_(),
      parameters_()
    {
    }

    /// copy constructor
    ModelDescription(const ModelDescription & source) :
      name_(source.name_),
      parameters_(source.parameters_)
    {
    }

    /// constructor provided for convenience
    ModelDescription(const BaseModel<D> * model) :
      name_(model->getName()),
      parameters_(model->getParameters())
    {
    }

    /// destructor
    virtual ~ModelDescription()
    {
    }

    /// assignment operator
    virtual ModelDescription & operator=(const ModelDescription & source)
    {
      if (&source == this) return *this;

      name_ = source.name_;
      parameters_ = source.parameters_;

      return *this;
    }

    /// creates model from the parameters defined in this class
    /// returns 0 if no description is set.
    BaseModel<D> * createModel()
    {
      if (name_ == "") return nullptr;

      BaseModel<D> * model = Factory<BaseModel<D> >::create(name_);
      model->setParameters(parameters_);
      return model;
    }

    /** Accessors */
    //@{
    /// Non-mutable access to model name
    const String & getName() const
    {
      return name_;
    }

    /// Mutable access to the model name
    String & getName()
    {
      return name_;
    }

    /// Set the model name
    void setName(const String & name)
    {
      name_ = name;
    }

    /// Non-mutable access to model parameters
    const Param & getParam() const
    {
      return parameters_;
    }

    /// Mutable access to the model parameters
    Param & getParam()
    {
      return parameters_;
    }

    /// Set the model parameters
    void setParam(const Param & param)
    {
      parameters_ = param;
    }

    /** @name Predicates */
    //@{
    virtual bool operator==(const ModelDescription & rhs) const
    {
      return (name_ == rhs.name_) && (parameters_ == rhs.parameters_);
    }

    virtual bool operator!=(const ModelDescription & rhs) const
    {
      return !(operator==(rhs));
    }

    //@}

protected:

    String name_;
    Param parameters_;
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_MODELDESCRIPTION_H
