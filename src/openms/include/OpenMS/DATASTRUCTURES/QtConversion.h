// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: Prachi Agrawal $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

// Only compile this header when building with GUI support
#ifdef WITH_GUI

#include <QString>
#include <QStringList>
#include <QByteArray>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
   * @brief Conversion utilities between OpenMS core types and Qt types
   * 
   * This header provides lightweight conversion functions to bridge the
   * Qt-free OpenMS core library with the Qt-dependent GUI components.
   * 
   * These functions are only available when building with -DWITH_GUI=ON
   * and should only be used at the boundary between core and GUI code.
   * 
   * @note Core library code should NEVER use these functions.
   * @note GUI library code should use these for conversion at boundaries only.
   * 
   * Example usage in GUI code:
   * @code
   * // Converting from core to GUI
   * std::string core_string = some_core_function();
   * QString qt_string = QtConversion::toQt(core_string);
   * 
   * // Converting from GUI to core
   * QString user_input = lineEdit->text();
   * std::string core_input = QtConversion::fromQt(user_input);
   * @endcode
   * 
   * @ingroup Datastructures
   */
  class OPENMS_DLLAPI QtConversion
  {
  public:
    
    ///@name QString conversions
    ///@{
    
    /**
     * @brief Convert std::string to QString
     * 
     * @param str Standard string to convert
     * @return QString UTF-8 encoded Qt string
     */
    static QString toQt(const std::string& str)
    {
      return QString::fromStdString(str);
    }
    
    /**
     * @brief Convert QString to std::string
     * 
     * @param qstr Qt string to convert
     * @return std::string UTF-8 encoded standard string
     */
    static std::string fromQt(const QString& qstr)
    {
      return qstr.toStdString();
    }
    
    ///@}
    
    ///@name QStringList conversions
    ///@{
    
    /**
     * @brief Convert std::vector<std::string> to QStringList
     * 
     * @param vec Vector of standard strings
     * @return QStringList Qt string list
     */
    static QStringList toQt(const std::vector<std::string>& vec)
    {
      QStringList result;
      result.reserve(static_cast<int>(vec.size()));
      for (const auto& str : vec)
      {
        result << QString::fromStdString(str);
      }
      return result;
    }
    
    /**
     * @brief Convert QStringList to std::vector<std::string>
     * 
     * @param qlist Qt string list
     * @return std::vector<std::string> Vector of standard strings
     */
    static std::vector<std::string> fromQt(const QStringList& qlist)
    {
      std::vector<std::string> result;
      result.reserve(static_cast<size_t>(qlist.size()));
      for (const auto& qstr : qlist)
      {
        result.push_back(qstr.toStdString());
      }
      return result;
    }
    
    ///@}
    
    ///@name QByteArray conversions
    ///@{
    
    /**
     * @brief Convert std::vector<char> to QByteArray
     * 
     * @param vec Vector of bytes
     * @return QByteArray Qt byte array
     */
    static QByteArray toQt(const std::vector<char>& vec)
    {
      return QByteArray(vec.data(), static_cast<int>(vec.size()));
    }
    
    /**
     * @brief Convert QByteArray to std::vector<char>
     * 
     * @param qarray Qt byte array
     * @return std::vector<char> Vector of bytes
     */
    static std::vector<char> fromQtToVector(const QByteArray& qarray)
    {
      return std::vector<char>(qarray.begin(), qarray.end());
    }
    
    /**
     * @brief Convert std::string to QByteArray
     * 
     * @param str Standard string (may contain binary data)
     * @return QByteArray Qt byte array
     */
    static QByteArray toQtBytes(const std::string& str)
    {
      return QByteArray(str.data(), static_cast<int>(str.size()));
    }
    
    /**
     * @brief Convert QByteArray to std::string
     * 
     * @param qarray Qt byte array
     * @return std::string Standard string (may contain binary data)
     */
    static std::string fromQtBytes(const QByteArray& qarray)
    {
      return std::string(qarray.constData(), static_cast<size_t>(qarray.size()));
    }
    
    ///@}
    
    ///@name Helper utilities
    ///@{
    
    /**
     * @brief Join a vector of strings with a delimiter (Qt-style)
     * 
     * @param vec Vector of strings to join
     * @param delimiter Delimiter to use between strings
     * @return QString Joined string
     */
    static QString joinToQt(const std::vector<std::string>& vec, const QString& delimiter)
    {
      return toQt(vec).join(delimiter);
    }
    
    /**
     * @brief Split a QString and convert to vector of strings
     * 
     * @param qstr String to split
     * @param delimiter Delimiter to split on
     * @return std::vector<std::string> Vector of split strings
     */
    static std::vector<std::string> splitFromQt(const QString& qstr, const QString& delimiter)
    {
      return fromQt(qstr.split(delimiter));
    }
    
    ///@}
    
  private:
    /// Private constructor - this is a utility class with static methods only
    QtConversion() = delete;
    ~QtConversion() = delete;
    QtConversion(const QtConversion&) = delete;
    QtConversion& operator=(const QtConversion&) = delete;
  };
  
} // namespace OpenMS

#endif // WITH_GUI
