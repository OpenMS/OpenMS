// Windows long-path probe -- INVESTIGATION ONLY, not part of the OpenMS build.
//
// Question: which of the file APIs OpenMS actually uses can handle a path longer
// than MAX_PATH (260) on Windows, and does a `longPathAware` manifest change the answer?
//
// Background. OpenMS 3.5 reached the filesystem through Qt. Qt's
// QFileSystemEntry::resolveNativeFilePath() (qfilesystementry.cpp:130) runs every path
// through QFSFileEnginePrivate::longFileName(), which makes it absolute and prepends
// "\\?\" -- so QFileInfo/QDir/QFile worked past 260 with no manifest at all.
// OpenMS 3.6 dropped Qt from the core (24 -> 1 files using QFile/QFileInfo/QDir) and
// uses std::filesystem and std::ifstream/ofstream. MS STL's __std_fs_open_handle()
// passes the path straight to CreateFileW with no prefixing, so those work past 260
// only if the process opted in (manifest + LongPathsEnabled registry value).
//
// Third-party libraries were never routed through Qt, so each decides for itself:
//   SQLite  os_win.c:2926-2937  prefixes "\\?\" itself when nChar > MAX_PATH   -> expected OK
//   Xerces  WindowsFileMgr.cpp:143  raw ::CreateFileW, no prefix               -> expected FAIL
//   Arrow   io_util.cc:1058/:1113   raw CreateFileW, no prefix                 -> expected FAIL
// The "CreateFileW (plain)" row below is exactly what Xerces and Arrow do, so it stands
// in for both without needing to link them.
//
// Build twice: once plain, once with longpath.manifest embedded. Compare.

#include <windows.h>
#include <io.h>       // _waccess_s, _wsopen_s, _close, _wunlink
#include <fcntl.h>    // _O_CREAT / _O_EXCL / _O_WRONLY
#include <share.h>    // _SH_DENYNO
#include <sys/stat.h> // _S_IREAD / _S_IWRITE
#include <cerrno>
#include <filesystem>
#include <fstream>
#include <cstdio>
#include <string>

#ifdef PROBE_WITH_SQLITE
#  include <sqlite3.h>
#endif

namespace fs = std::filesystem;

static int g_pass = 0, g_fail = 0;

static void row(const char* api, const char* models, bool ok, unsigned long err)
{
    printf("  %-28s | %-34s | %-6s | %lu\n", api, models, ok ? "OK" : "FAIL", ok ? 0UL : err);
    ok ? ++g_pass : ++g_fail;
}

int main()
{
    printf("MAX_PATH = %d\n", MAX_PATH);

    // Is this process long-path aware? Only observable indirectly; report the registry value.
    {
        HKEY k; DWORD v = 0, sz = sizeof(v), type = 0;
        if (RegOpenKeyExW(HKEY_LOCAL_MACHINE,
                L"SYSTEM\\CurrentControlSet\\Control\\FileSystem", 0, KEY_READ, &k) == ERROR_SUCCESS) {
            RegQueryValueExW(k, L"LongPathsEnabled", nullptr, &type, (LPBYTE)&v, &sz);
            RegCloseKey(k);
        }
        printf("HKLM\\...\\FileSystem\\LongPathsEnabled = %lu\n", v);
    }
#ifdef PROBE_MANIFEST
    printf("manifest: longPathAware = true (embedded)\n");
#else
    printf("manifest: none\n");
#endif

    // Build a directory tree past MAX_PATH. Only the \\?\ form can create it.
    wchar_t cwd[32768]; GetCurrentDirectoryW(32768, cwd);
    std::wstring base = std::wstring(cwd) + L"\\lp";
    CreateDirectoryW((L"\\\\?\\" + base).c_str(), nullptr);
    std::wstring deep = base;
    const std::wstring comp(40, L'd');
    for (int i = 0; i < 12; ++i) {
        std::wstring next = deep + L"\\" + comp;
        if (!CreateDirectoryW((L"\\\\?\\" + next).c_str(), nullptr)
            && GetLastError() != ERROR_ALREADY_EXISTS) break;
        deep = next;
    }
    const std::wstring wtarget = deep + L"\\out.dat";
    std::string ntarget; ntarget.reserve(wtarget.size());         // ASCII by construction
    for (wchar_t c : wtarget) ntarget.push_back(static_cast<char>(c));
    printf("target path length = %zu\n\n", wtarget.size());

    printf("  %-28s | %-34s | %-6s | %s\n", "API", "stands in for", "result", "GetLastError/errno");
    printf("  %s\n", std::string(28, '-').append("-+-").append(34, '-')
                     .append("-+-").append(6, '-').append("-+-").append(18, '-').c_str());

    // --- what OpenMS 3.6 uses ------------------------------------------------
    { std::ofstream os(ntarget.c_str());
      row("std::ofstream(const char*)", "OpenMS raw streams (255 sites)", os.is_open() && os.good(), (unsigned long)errno); }

    { std::ofstream os{fs::path(wtarget)};   // braces: (fs::path(x)) would be a function declaration
      row("std::ofstream(fs::path)", "OpenMS via std::filesystem", os.is_open() && os.good(), (unsigned long)errno); }

    { std::error_code ec; bool ok = fs::create_directories(fs::path(deep) / L"sub", ec) || !ec;
      row("fs::create_directories", "File::TempDir / mkdir paths", ok, (unsigned long)ec.value()); }

    { errno = 0; bool ok = (_waccess_s(wtarget.c_str(), 0) == 0);
      row("_waccess_s(existence)", "File::exists / readable / writable", ok, (unsigned long)errno); }

    // --- what Xerces and Arrow do (raw CreateFileW, no prefix) ---------------
    { HANDLE h = CreateFileW(wtarget.c_str(), GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
      bool ok = h != INVALID_HANDLE_VALUE; DWORD e = ok ? 0 : GetLastError(); if (ok) CloseHandle(h);
      row("CreateFileW (plain)", "Xerces mzML / Arrow Parquet", ok, e); }

    // --- what Qt did in 3.5, and what SQLite does for itself ------------------
    { std::wstring p = L"\\\\?\\" + wtarget;
      HANDLE h = CreateFileW(p.c_str(), GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
      bool ok = h != INVALID_HANDLE_VALUE; DWORD e = ok ? 0 : GetLastError(); if (ok) CloseHandle(h);
      row("CreateFileW (\\\\?\\ prefix)", "OpenMS 3.5 via Qt longFileName()", ok, e); }

    // --- would prefixing at a chokepoint restore 3.5 parity? -----------------
    // File::readable/writable use the NARROW _access_s/_sopen_s today. Parity would mean
    // switching to the wide forms fed a \\?\ path. These rows say whether the CRT wide
    // functions actually honour the prefix -- the load-bearing unknown for that plan.
    { const std::wstring pp = L"\\\\?\\" + wtarget;
      errno = 0; bool ok = (_waccess_s(pp.c_str(), 0) == 0);
      row("_waccess_s (\\\\?\\ prefix)", "parity fix: File::readable/writable", ok, (unsigned long)errno); }

    { const std::wstring pp = L"\\\\?\\" + wtarget + L".excl";
      int fd = -1; errno = 0;
      const int e = _wsopen_s(&fd, pp.c_str(), _O_CREAT | _O_EXCL | _O_WRONLY, _SH_DENYNO, _S_IREAD | _S_IWRITE);
      bool ok = (e == 0 && fd >= 0);
      if (fd >= 0) _close(fd);
      _wunlink(pp.c_str());   // must delete: both binaries share a cwd, and _O_EXCL
                              // would otherwise report EEXIST on the second run
      row("_wsopen_s (\\\\?\\ prefix)", "parity fix: writable() probe", ok, (unsigned long)(ok ? 0 : errno)); }

    { std::ofstream os{fs::path(L"\\\\?\\" + wtarget)};
      row("std::ofstream (\\\\?\\ path)", "parity fix: raw streams", os.is_open() && os.good(), (unsigned long)errno); }

#ifdef PROBE_WITH_SQLITE
    { std::string db = ntarget + ".sqlite";
      sqlite3* h = nullptr;
      int rc = sqlite3_open(db.c_str(), &h);
      bool ok = (rc == SQLITE_OK);
      if (ok) { char* err = nullptr;
                ok = sqlite3_exec(h, "create table t(x); insert into t values(1);", nullptr, nullptr, &err) == SQLITE_OK;
                if (err) sqlite3_free(err); }
      if (h) sqlite3_close(h);
      row("sqlite3_open + write", "OpenMS SQLite (72 sites)", ok, (unsigned long)rc); }
#endif

    printf("\n  %d passed, %d failed\n", g_pass, g_fail);
    return 0;  // never fail the job; this is a measurement, not an assertion
}
