\# vcpkg Integration Verification (Windows)



\## Environment



\* OS: Windows 11

\* Compiler: MSVC (Visual Studio 2022 Build Tools)

\* CMake: 4.3.1

\* vcpkg: latest (GitHub)



\## Steps Performed



1\. Installed dependencies using `vcpkg.json`

2\. Configured OpenMS using CMake with vcpkg toolchain

3\. Built OpenMS successfully





\## Issues Encountered



\### 1. Toolchain Not Applied (Critical)



Using relative path:



\-DCMAKE\_TOOLCHAIN\_FILE=../vcpkg/scripts/buildsystems/vcpkg.cmake



Result:



\* CMake ignored toolchain

\* Dependencies like XercesC not found



Fix:



\-DCMAKE\_TOOLCHAIN\_FILE="C:\\absolute\\path\\to\\vcpkg.cmake"



\### 2. Silent Failure of CMake



CMake does not error when toolchain is ignored, leading to confusing dependency errors.



\### 3. Missing dumpbin.exe (Windows-specific)



Error:



Could not find 'dumpbin.exe'



Cause:



\* MSVC tools not in PATH



Fix:



\* Use \*\*"x64 Native Tools Command Prompt for VS 2022"\*\*





\## Final Result



\* All dependencies installed successfully

\* OpenMS builds successfully using vcpkg

\* Verified working configuration on Windows



\## Suggestions



\* Recommend absolute toolchain path in documentation

\* Add warning when toolchain is ignored

\* Document requirement of VS Developer Command Prompt on Windows



