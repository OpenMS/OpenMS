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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MetaInfoInterface.h>

using namespace std;

namespace OpenMS
{

  MetaInfoInterface::MetaInfoInterface() :
    meta_(nullptr)
  {

  }

  MetaInfoInterface::MetaInfoInterface(const MetaInfoInterface & rhs)
  {
    if (rhs.meta_ != nullptr)
    {
      meta_ = new MetaInfo(*(rhs.meta_));
    }
    else
    {
      meta_ = nullptr;
    }
  }

  MetaInfoInterface::~MetaInfoInterface()
  {
    delete(meta_);
  }

  MetaInfoInterface & MetaInfoInterface::operator=(const MetaInfoInterface & rhs)
  {
    if (this == &rhs)
      return *this;

//      std::cout << meta_ << std::endl;
//      std::cout << rhs.meta_ << std::endl;
//      std::cout << " " << std::endl;

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

  bool MetaInfoInterface::operator==(const MetaInfoInterface & rhs) const
  {
    if (rhs.meta_ == nullptr && meta_ == nullptr)
    {
      return true;
    }
    else if (rhs.meta_ == nullptr && meta_ != nullptr)
    {
      if (meta_->empty())
        return true;

      return false;
    }
    else if (rhs.meta_ != nullptr && meta_ == nullptr)
    {
      if (rhs.meta_->empty())
        return true;

      return false;
    }
    return *meta_ == *(rhs.meta_);
  }

  bool MetaInfoInterface::operator!=(const MetaInfoInterface & rhs) const
  {
    return !(operator==(rhs));
  }

  const DataValue & MetaInfoInterface::getMetaValue(const String & name) const
  {
    if (meta_ == nullptr)
    {
      return DataValue::EMPTY;
    }
    return meta_->getValue(name);
  }

  const DataValue & MetaInfoInterface::getMetaValue(UInt index) const
  {
    if (meta_ == nullptr)
    {
      return DataValue::EMPTY;
    }
    return meta_->getValue(index);
  }

  bool MetaInfoInterface::metaValueExists(const String & name) const
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

  void MetaInfoInterface::setMetaValue(const String & name, const DataValue & value)
  {
    createIfNotExists_();
    meta_->setValue(name, value);
  }

  void MetaInfoInterface::setMetaValue(UInt index, const DataValue & value)
  {
    createIfNotExists_();
    meta_->setValue(index, value);
  }

  MetaInfoRegistry & MetaInfoInterface::metaRegistry()
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

  void MetaInfoInterface::getKeys(std::vector<String> & keys) const
  {
    if (meta_ != nullptr)
    {
      meta_->getKeys(keys);
    }
  }

  void MetaInfoInterface::getKeys(std::vector<UInt> & keys) const
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

  void MetaInfoInterface::removeMetaValue(const String & name)
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

} //namespace
