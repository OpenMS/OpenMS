// src/openms/thirdparty/percolator/Logger.h  (OpenMS shim — original contents
// replaced by patch 06-logger-bridge.patch)
#ifndef LOGGER_H_
#define LOGGER_H_

#include <OpenMS/CONCEPT/LogStream.h>
#include <sstream>


namespace OpenMS { namespace Internal { namespace Percolator {
enum class LogLevel { INFO, DEBUG, WARN, ERROR };

class Log
{
public:
  Log(LogLevel lvl) : level_(lvl) {}
  ~Log()
  {
    const std::string msg = buffer_.str();
    switch (level_)
    {
      case LogLevel::INFO:  OPENMS_LOG_INFO  << "[percolator] " << msg << '\n'; break;
      case LogLevel::DEBUG: OPENMS_LOG_DEBUG << "[percolator] " << msg << '\n'; break;
      case LogLevel::WARN:  OPENMS_LOG_WARN  << "[percolator] " << msg << '\n'; break;
      case LogLevel::ERROR: OPENMS_LOG_ERROR << "[percolator] " << msg << '\n'; break;
    }
  }
  template<typename T> Log& operator<<(const T& v) { buffer_ << v; return *this; }
private:
  LogLevel level_;
  std::ostringstream buffer_;
};

// Legacy macros used by upstream Percolator code. Re-route to OpenMS log streams.
#define LOG_INFO  ::OpenMS::Internal::Percolator::Log(::OpenMS::Internal::Percolator::LogLevel::INFO)
#define LOG_DEBUG ::OpenMS::Internal::Percolator::Log(::OpenMS::Internal::Percolator::LogLevel::DEBUG)
#define LOG_WARN  ::OpenMS::Internal::Percolator::Log(::OpenMS::Internal::Percolator::LogLevel::WARN)
#define LOG_ERROR ::OpenMS::Internal::Percolator::Log(::OpenMS::Internal::Percolator::LogLevel::ERROR)

}}}  // namespace OpenMS::Internal::Percolator
#endif
