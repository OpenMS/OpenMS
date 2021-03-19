//
// Created by david on 18.03.21.
//

#pragma once

#include <OpenMS/OpenMSConfig.h>

#include <future>
#include <vector>

namespace OpenMS
{
    class Param;
    class String;

    class OPENMS_DLLAPI TVToolDiscovery
    {
    private:
      static Param getParamFromIni_(const String& tool_name);

    public:
      static std::vector<std::future<Param>> results_;

      static void findTools();

      static std::vector<Param> getToolParams();
    };
}