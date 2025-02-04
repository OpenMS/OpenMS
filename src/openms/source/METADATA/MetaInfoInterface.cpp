// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <OpenMS/METADATA/MetaInfo.h>

using namespace std;

namespace OpenMS
{

  MetaInfoInterface::MetaInfoInterface() :
    meta_(nullptr)
  {
  }

  /// Copy constructor
  MetaInfoInterface::MetaInfoInterface(const MetaInfoInterface& rhs) :
    meta_(nullptr)
  {
    if (rhs.meta_ != nullptr)
    {
      meta_ = new MetaInfo(*(rhs.meta_));
    }
  }

  /// Move constructor
  MetaInfoInterface::MetaInfoInterface(MetaInfoInterface&& rhs) noexcept :
    meta_(std::move(rhs.meta_))
  {
    // take ownership
    rhs.meta_ = nullptr;
  }

  MetaInfoInterface::~MetaInfoInterface()
  {
    delete(meta_);
  }

  MetaInfoInterface& MetaInfoInterface::operator=(const MetaInfoInterface& rhs)
  {
    if (this == &rhs)
    {
      return *this;
    }

    if (rhs.meta_ != nullptr && meta_ != nullptr)
    {
      *meta_ = *(rhs.meta_);
    }
    else if (rhs.meta_ == nullptr && meta_ != nullptr)
    {
      delete(meta_);
      meta_ = nullptr;
    }
    else if (rhs.meta_ != nullptr && meta_ == nullptr)
    {
      meta_ = new MetaInfo(*(rhs.meta_));
    }

    return *this;
  }

  // Move assignment
  MetaInfoInterface& MetaInfoInterface::operator=(MetaInfoInterface&& rhs) noexcept
  {
    if (this == &rhs)
    {
      return *this;
    }

    // free memory and assign rhs memory
    delete(meta_);
    meta_ = rhs.meta_;
    rhs.meta_ = nullptr;

    return *this;
  }

  void MetaInfoInterface::swap(MetaInfoInterface& rhs)
  {
    std::swap(meta_, rhs.meta_);
  //   MetaInfo* temp = meta_;
  //   meta_ = rhs.meta_;
  //   rhs.meta_ = temp;
  }

  bool MetaInfoInterface::operator==(const MetaInfoInterface& rhs) const
  {
    if (rhs.meta_ == nullptr && meta_ == nullptr)
    {
      return true;
    }
    else if (rhs.meta_ == nullptr && meta_ != nullptr)
    {
      if (meta_->empty())
      {
        return true;
      }
      return false;
    }
    else if (rhs.meta_ != nullptr && meta_ == nullptr)
    {
      if (rhs.meta_->empty())
      {
        return true;
      }
      return false;
    }
    return *meta_ == *(rhs.meta_);
  }

  bool MetaInfoInterface::operator!=(const MetaInfoInterface& rhs) const
  {
    return !(operator==(rhs));
  }

  const DataValue& MetaInfoInterface::getMetaValue(const String& name) const
  {
    if (meta_ == nullptr)
    {
      return DataValue::EMPTY;
    }
    return meta_->getValue(name, DataValue::EMPTY);
  }

  DataValue MetaInfoInterface::getMetaValue(const String& name, const DataValue& default_value) const
  {
    if (meta_ == nullptr)
    {
      return default_value;
    }
    return meta_->getValue(name, default_value);
  }

  const DataValue& MetaInfoInterface::getMetaValue(UInt index) const
  {
    if (meta_ == nullptr)
    {
      return DataValue::EMPTY;
    }
    return meta_->getValue(index, DataValue::EMPTY);
  }

  DataValue MetaInfoInterface::getMetaValue(UInt index, const DataValue& default_value) const
  {
    if (meta_ == nullptr)
    {
      return default_value;
    }
    return meta_->getValue(index, default_value);
  }

  bool MetaInfoInterface::metaValueExists(const String& name) const
  {
    if (meta_ == nullptr)
    {
      return false;
    }
    return meta_->exists(name);
  }

  bool MetaInfoInterface::metaValueExists(UInt index) const
  {
    if (meta_ == nullptr)
    {
      return false;
    }
    return meta_->exists(index);
  }

  void MetaInfoInterface::setMetaValue(const String& name, const DataValue& value)
  {
    createIfNotExists_();
    meta_->setValue(name, value);
  }

  void MetaInfoInterface::setMetaValue(UInt index, const DataValue& value)
  {
    createIfNotExists_();
    meta_->setValue(index, value);
  }

  MetaInfoRegistry& MetaInfoInterface::metaRegistry()
  {
    return MetaInfo::registry();
  }

  void MetaInfoInterface::createIfNotExists_()
  {
    if (meta_ == nullptr)
    {
      meta_ = new MetaInfo();
    }
  }

  void MetaInfoInterface::getKeys(std::vector<String>& keys) const
  {
    if (meta_ != nullptr)
    {
      meta_->getKeys(keys);
    }
  }

  void MetaInfoInterface::getKeys(std::vector<UInt>& keys) const
  {
    if (meta_ != nullptr)
    {
      meta_->getKeys(keys);
    }
  }

  bool MetaInfoInterface::isMetaEmpty() const
  {
    if (meta_ == nullptr)
    {
      return true;
    }
    return meta_->empty();
  }

  void MetaInfoInterface::clearMetaInfo()
  {
    delete meta_;
    meta_ = nullptr;
  }

  void MetaInfoInterface::removeMetaValue(const String& name)
  {
    if (meta_ != nullptr)
    {
      meta_->removeValue(name);
    }
  }

  void MetaInfoInterface::removeMetaValue(UInt index)
  {
    if (meta_ != nullptr)
    {
      meta_->removeValue(index);
    }
  }

  //TODO get a MetaValue list to copy only those that have been set
  void MetaInfoInterface::addMetaValues(const MetaInfoInterface& from)
  {
    std::vector<String> keys;
    from.getKeys(keys);
    for (String& key : keys)
    {
      this->setMetaValue(key, from.getMetaValue(key));
    }
  }

} //namespace
