/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/


#include "Logger.h"

Logger::Logger()
{
  file_log = false;
}

Logger::Logger(const char* file)
{
  attach_file(file);
}

void Logger::attach_file(const char* file)
{
  std_log.flush();
  std_log.close();
  std_log.open(file);
  if (!(std_log.is_open())) {
       std::cerr << "[]: Couldn't open file \"" << file << "\" for logging.\n";
       this->~Logger(); /* Destroy the object */
  }
  file_log = true;
}

Logger::~Logger()
{
  std::cerr.flush();
  std_log.close();
}
