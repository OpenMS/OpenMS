// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FEATUREFINDER/BaseModel.h>

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
    ModelDescription(const BaseModel * model) :
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
