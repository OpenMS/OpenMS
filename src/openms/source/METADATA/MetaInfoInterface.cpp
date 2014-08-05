// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <QtCore/QList>
#include <QtCore/QStringList>
#include <QtCore/QVariant>

using namespace std;

namespace OpenMS
{

  MetaInfoInterface::MetaInfoInterface() :
    meta_(0)
  {

  }

  MetaInfoInterface::MetaInfoInterface(const MetaInfoInterface & rhs)
  {
    if (rhs.meta_ != 0)
    {
      meta_ = new MetaInfo(*(rhs.meta_));
    }
    else
    {
      meta_ = 0;
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

    if (rhs.meta_ != 0 && meta_ != 0)
    {
      *meta_ = *(rhs.meta_);
    }
    else if (rhs.meta_ == 0 && meta_ != 0)
    {
      delete(meta_);
      meta_ = 0;
    }
    else if (rhs.meta_ != 0 && meta_ == 0)
    {
      meta_ = new MetaInfo(*(rhs.meta_));
    }

    return *this;
  }

  bool MetaInfoInterface::operator==(const MetaInfoInterface & rhs) const
  {
    if (rhs.meta_ == 0 && meta_ == 0)
    {
      return true;
    }
    else if (rhs.meta_ == 0 && meta_ != 0)
    {
      if (meta_->empty())
        return true;

      return false;
    }
    else if (rhs.meta_ != 0 && meta_ == 0)
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
    if (meta_ == 0)
    {
      return DataValue::EMPTY;
    }
    return meta_->getValue(name);
  }

  const DataValue & MetaInfoInterface::getMetaValue(UInt index) const
  {
    if (meta_ == 0)
    {
      return DataValue::EMPTY;
    }
    return meta_->getValue(index);
  }

  bool MetaInfoInterface::metaValueExists(const String & name) const
  {
    if (meta_ == 0)
    {
      return false;
    }
    return meta_->exists(name);
  }

  bool MetaInfoInterface::metaValueExists(UInt index) const
  {
    if (meta_ == 0)
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
    if (meta_ == 0)
    {
      meta_ = new MetaInfo();
    }
  }

  void MetaInfoInterface::getKeys(std::vector<String> & keys) const
  {
    if (meta_ != 0)
    {
      meta_->getKeys(keys);
    }
  }

  void MetaInfoInterface::getKeys(std::vector<UInt> & keys) const
  {
    if (meta_ != 0)
    {
      meta_->getKeys(keys);
    }
  }

  bool MetaInfoInterface::isMetaEmpty() const
  {
    if (meta_ == 0)
    {
      return true;
    }
    return meta_->empty();
  }

  void MetaInfoInterface::clearMetaInfo()
  {
    delete meta_;
    meta_ = 0;
  }

  void MetaInfoInterface::removeMetaValue(const String & name)
  {
    if (meta_ != 0)
    {
      meta_->removeValue(name);
    }
  }

  void MetaInfoInterface::removeMetaValue(UInt index)
  {
    if (meta_ != 0)
    {
      meta_->removeValue(index);
    }
  }

  /// In-stream operator to serialize MetaInfoInterface instance
  OPENMS_DLLAPI QDataStream& operator>>(QDataStream& in, MetaInfoInterface& metaInfoInterface)
  {
    QList<QVariant> keysList;
    QList<QVariant> typeList;
    QList<QVariant> valueList;
    QList<QVariant> unitsList;

    QString temp;
    QStringList tempList;
    StringList sList;
    IntList iList;
    DoubleList dList;

    in >> keysList >> typeList >> valueList >> unitsList;

    for (int i = 0; i < keysList.size(); i++)
    {
      String key = keysList[i].toString();

      DataValue dv;
      switch(typeList[i].toInt())
      {
      case DataValue::STRING_VALUE:
          dv = valueList[i].toString();
          break;

      case DataValue::INT_VALUE:
          dv = valueList[i].toString().toInt();
          break;

      case DataValue::DOUBLE_VALUE:
          dv = valueList[i].toString().toDouble();
          break;

      case DataValue::STRING_LIST:
          temp = valueList[i].toString().mid(1); // remove left quote
          temp.chop(1); //remove right quote
          tempList = temp.split(',');
          for (int j = 0; j < tempList.size(); j++)
          {
              sList.push_back(String(tempList[j].trimmed()));
          }
          dv = sList;
          break;

      case DataValue::INT_LIST:
          temp = valueList[i].toString().mid(1); // remove left quote
          temp.chop(1); //remove right quote
          tempList = temp.split(',');

          for (int j = 0; j < tempList.size(); j++)
          {
              iList.push_back(tempList[j].trimmed().toInt());
          }
          dv = iList;
          break;

      case DataValue::DOUBLE_LIST:
          temp = valueList[i].toString().mid(1); // remove left quote
          temp.chop(1); //remove right quote
          tempList = temp.split(',');
          for (int j = 0; j < tempList.size(); j++)
          {
              dList.push_back(tempList[j].trimmed().toDouble());
          }
          dv = dList;
          break;
      }

      if(unitsList[i].isValid())
      {
          dv.setUnit(unitsList[i].toString());
      }

      metaInfoInterface.setMetaValue(key, dv);

    }

    return in;

  }

  /// Out-stream operator to serialize MetaInfoInterface instance
  OPENMS_DLLAPI QDataStream& operator<<(QDataStream& out, const MetaInfoInterface& metaInfoInterface)
  {
    vector<String> keys;
    metaInfoInterface.getKeys(keys);

    QList<QVariant> keysList;
    QList<QVariant> typeList;
    QList<QVariant> valueList;
    QList<QVariant> unitsList;

    vector<String>::const_iterator iter;
    for (iter = keys.begin(); iter != keys.end(); ++iter)
    {
      String key = *iter;

      keysList.append(key.toQString());

      DataValue dv = metaInfoInterface.getMetaValue(key);
      QVariant type = QVariant((int) dv.valueType());
      typeList.append(type);

      valueList.append(dv.toQString());

      if(dv.hasUnit())
      {
        unitsList.append(dv.getUnit().toQString());
      }
      else
      {
        unitsList.append(QVariant());
      }

    }

    out << keysList << typeList << valueList << unitsList;

    return out;

  }


} //namespace
