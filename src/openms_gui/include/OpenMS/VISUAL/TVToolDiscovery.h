//
// Created by david on 18.03.21.
//

#pragma once

#include <OpenMS/OpenMSConfig.h>

#include <future>
#include <vector>
#include <unordered_map>

namespace OpenMS
{
    class Param;
    class String;

    class OPENMS_DLLAPI TVToolDiscovery
    {
    private:
      static Param getParamFromIni_(const String& tool_name);

      static std::unordered_map<std::string, std::future<Param>> future_results_;

      static bool ready_;

      static std::unordered_map<std::string, Param> params_;

    public:
      static void findTools();

      static std::unordered_map<std::string, Param>& getToolParams();
    };
}