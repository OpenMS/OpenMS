// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/ExternalProcess.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <cstdlib>
#include <utility>
#include <sstream>
#include <thread>
#include <chrono>
#include <iostream>

#ifdef OPENMS_WINDOWSPLATFORM
  #include <windows.h>
  #include <io.h>
  #include <fcntl.h>
#else
  #include <unistd.h>
  #include <sys/wait.h>
  #include <sys/types.h>
  #include <fcntl.h>
  #include <cerrno>
  #include <cstring>
  #include <signal.h>
#endif

namespace OpenMS
{

  /// default Ctor; callbacks for stdout/stderr are empty
  ExternalProcess::ExternalProcess()
    : ExternalProcess([](const String& /*out*/) {}, [](const String& /*out*/) {})
  {
  }

  ExternalProcess::ExternalProcess(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr)
    : callbackStdOut_(std::move(callbackStdOut)),
      callbackStdErr_(std::move(callbackStdErr))
  {
  }

  ExternalProcess::~ExternalProcess()
  {
  }

  /// re-wire the callbacks used using run()
  void ExternalProcess::setCallbacks(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr)
  {
    callbackStdOut_ = std::move(callbackStdOut);
    callbackStdErr_ = std::move(callbackStdErr);
  }


  ExternalProcess::RETURNSTATE ExternalProcess::run(const String& exe, const std::vector<String>& args, const String& working_dir, const bool verbose, IO_MODE io_mode, const std::map<String, String>& env)
  {
    String error_msg;
    return run(exe, args, working_dir, verbose, error_msg, io_mode, env);
  }

  ExternalProcess::RETURNSTATE ExternalProcess::run(const String& exe, const std::vector<String>& args, const String& working_dir, const bool verbose, String& error_msg, IO_MODE /*io_mode*/, const std::map<String, String>& env)
  {
    error_msg.clear();

    // Build command line
    std::ostringstream cmd_stream;
    cmd_stream << exe;
    for (const auto& arg : args)
    {
      cmd_stream << " " << arg;
    }
    String full_cmd = cmd_stream.str();

    if (verbose)
    {
      callbackStdOut_(String("Running: ") + full_cmd + '\n');
    }

#ifdef OPENMS_WINDOWSPLATFORM
    // Windows implementation using CreateProcess
    SECURITY_ATTRIBUTES saAttr;
    saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
    saAttr.bInheritHandle = TRUE;
    saAttr.lpSecurityDescriptor = NULL;

    HANDLE hChildStdOutRd = NULL, hChildStdOutWr = NULL;
    HANDLE hChildStdErrRd = NULL, hChildStdErrWr = NULL;

    // Create pipes for stdout and stderr
    if (!CreatePipe(&hChildStdOutRd, &hChildStdOutWr, &saAttr, 0))
    {
      error_msg = "Failed to create stdout pipe";
      if (verbose) callbackStdErr_(error_msg + '\n');
      return RETURNSTATE::FAILED_TO_START;
    }
    SetHandleInformation(hChildStdOutRd, HANDLE_FLAG_INHERIT, 0);

    if (!CreatePipe(&hChildStdErrRd, &hChildStdErrWr, &saAttr, 0))
    {
      error_msg = "Failed to create stderr pipe";
      if (verbose) callbackStdErr_(error_msg + '\n');
      CloseHandle(hChildStdOutRd);
      CloseHandle(hChildStdOutWr);
      return RETURNSTATE::FAILED_TO_START;
    }
    SetHandleInformation(hChildStdErrRd, HANDLE_FLAG_INHERIT, 0);

    PROCESS_INFORMATION piProcInfo;
    STARTUPINFOA siStartInfo;
    ZeroMemory(&piProcInfo, sizeof(PROCESS_INFORMATION));
    ZeroMemory(&siStartInfo, sizeof(STARTUPINFO));

    siStartInfo.cb = sizeof(STARTUPINFO);
    siStartInfo.hStdError = hChildStdErrWr;
    siStartInfo.hStdOutput = hChildStdOutWr;
    siStartInfo.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    siStartInfo.dwFlags |= STARTF_USESTDHANDLES;

    // Set environment variables
    char* envBlock = nullptr;
    if (!env.empty())
    {
      std::ostringstream envStream;
      for (const auto& kv : env)
      {
        envStream << kv.first << "=" << kv.second << '\0';
      }
      envStream << '\0';
      String envStr = envStream.str();
      envBlock = new char[envStr.size()];
      memcpy(envBlock, envStr.c_str(), envStr.size());
    }

    // Set working directory
    const char* workdir_cstr = working_dir.empty() ? nullptr : working_dir.c_str();

    // CreateProcessA may modify the command line buffer, so make a mutable copy
    std::vector<char> cmd_buffer(full_cmd.begin(), full_cmd.end());
    cmd_buffer.push_back('\0'); // Ensure null termination

    BOOL bSuccess = CreateProcessA(
      NULL,
      cmd_buffer.data(),
      NULL, NULL,
      TRUE,
      0,
      envBlock,
      workdir_cstr,
      &siStartInfo,
      &piProcInfo
    );

    delete[] envBlock;

    if (!bSuccess)
    {
      error_msg = "Process '" + exe + "' failed to start. Does it exist? Is it executable?";
      if (verbose) callbackStdErr_(error_msg + '\n');
      CloseHandle(hChildStdOutRd);
      CloseHandle(hChildStdOutWr);
      CloseHandle(hChildStdErrRd);
      CloseHandle(hChildStdErrWr);
      return RETURNSTATE::FAILED_TO_START;
    }

    CloseHandle(hChildStdOutWr);
    CloseHandle(hChildStdErrWr);

    // Read output
    const int BUFSIZE = 4096;
    char buffer[BUFSIZE];
    DWORD dwRead;

    while (true)
    {
      DWORD waitResult = WaitForSingleObject(piProcInfo.hProcess, 50);

      // Read stdout
      DWORD bytesAvail = 0;
      PeekNamedPipe(hChildStdOutRd, NULL, 0, NULL, &bytesAvail, NULL);
      if (bytesAvail > 0)
      {
        if (ReadFile(hChildStdOutRd, buffer, BUFSIZE - 1, &dwRead, NULL) && dwRead > 0)
        {
          buffer[dwRead] = '\0';
          callbackStdOut_(String(buffer));
        }
      }

      // Read stderr
      bytesAvail = 0;
      PeekNamedPipe(hChildStdErrRd, NULL, 0, NULL, &bytesAvail, NULL);
      if (bytesAvail > 0)
      {
        if (ReadFile(hChildStdErrRd, buffer, BUFSIZE - 1, &dwRead, NULL) && dwRead > 0)
        {
          buffer[dwRead] = '\0';
          callbackStdErr_(String(buffer));
        }
      }

      if (waitResult == WAIT_OBJECT_0) break;
    }

    // Final read to get any remaining output
    while (ReadFile(hChildStdOutRd, buffer, BUFSIZE - 1, &dwRead, NULL) && dwRead > 0)
    {
      buffer[dwRead] = '\0';
      callbackStdOut_(String(buffer));
    }
    while (ReadFile(hChildStdErrRd, buffer, BUFSIZE - 1, &dwRead, NULL) && dwRead > 0)
    {
      buffer[dwRead] = '\0';
      callbackStdErr_(String(buffer));
    }

    DWORD exitCode;
    GetExitCodeProcess(piProcInfo.hProcess, &exitCode);

    CloseHandle(piProcInfo.hProcess);
    CloseHandle(piProcInfo.hThread);
    CloseHandle(hChildStdOutRd);
    CloseHandle(hChildStdErrRd);

    if (exitCode != 0)
    {
      error_msg = "Process '" + exe + "' did not finish successfully (exit code: " + std::to_string(exitCode) + "). Please check the log.";
      if (verbose) callbackStdErr_(error_msg + '\n');
      return RETURNSTATE::NONZERO_EXIT;
    }

#else
    // Unix/Linux/macOS implementation using fork/exec
    int stdout_pipe[2], stderr_pipe[2];

    if (pipe(stdout_pipe) == -1 || pipe(stderr_pipe) == -1)
    {
      error_msg = "Failed to create pipes";
      if (verbose) callbackStdErr_(error_msg + '\n');
      return RETURNSTATE::FAILED_TO_START;
    }

    pid_t pid = fork();

    if (pid == -1)
    {
      error_msg = "Failed to fork process";
      if (verbose) callbackStdErr_(error_msg + '\n');
      close(stdout_pipe[0]);
      close(stdout_pipe[1]);
      close(stderr_pipe[0]);
      close(stderr_pipe[1]);
      return RETURNSTATE::FAILED_TO_START;
    }

    if (pid == 0)
    {
      // Child process
      close(stdout_pipe[0]);
      close(stderr_pipe[0]);

      dup2(stdout_pipe[1], STDOUT_FILENO);
      dup2(stderr_pipe[1], STDERR_FILENO);

      close(stdout_pipe[1]);
      close(stderr_pipe[1]);

      // Change working directory if specified
      if (!working_dir.empty())
      {
        if (chdir(working_dir.c_str()) != 0)
        {
          std::cerr << "Failed to change directory to: " << working_dir << std::endl;
          _exit(127);
        }
      }

      // Set environment variables in child only (not parent)
      for (const auto& kv : env)
      {
        setenv(kv.first.c_str(), kv.second.c_str(), 1);
      }

      // Build argument array
      std::vector<char*> argv_ptrs;
      argv_ptrs.push_back(const_cast<char*>(exe.c_str()));
      for (const auto& arg : args)
      {
        argv_ptrs.push_back(const_cast<char*>(arg.c_str()));
      }
      argv_ptrs.push_back(nullptr);

      execvp(exe.c_str(), argv_ptrs.data());

      // If we get here, exec failed
      std::cerr << "Failed to execute: " << exe << " - " << strerror(errno) << std::endl;
      _exit(127);
    }

    // Parent process
    close(stdout_pipe[1]);
    close(stderr_pipe[1]);

    // Set pipes to non-blocking
    fcntl(stdout_pipe[0], F_SETFL, O_NONBLOCK);
    fcntl(stderr_pipe[0], F_SETFL, O_NONBLOCK);

    // Read output
    const int BUFSIZE = 4096;
    char buffer[BUFSIZE];
    bool process_running = true;
    int status = 0;

    while (process_running)
    {
      // Check if process is still running
      pid_t result = waitpid(pid, &status, WNOHANG);

      if (result == pid)
      {
        process_running = false;
      }

      // Read stdout
      ssize_t bytes_read;
      while ((bytes_read = read(stdout_pipe[0], buffer, BUFSIZE - 1)) > 0)
      {
        buffer[bytes_read] = '\0';
        callbackStdOut_(String(buffer));
      }

      // Read stderr
      while ((bytes_read = read(stderr_pipe[0], buffer, BUFSIZE - 1)) > 0)
      {
        buffer[bytes_read] = '\0';
        callbackStdErr_(String(buffer));
      }

      if (process_running)
      {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
      }
    }

    // Final read to get any remaining output
    ssize_t bytes_read;
    while ((bytes_read = read(stdout_pipe[0], buffer, BUFSIZE - 1)) > 0)
    {
      buffer[bytes_read] = '\0';
      callbackStdOut_(String(buffer));
    }
    while ((bytes_read = read(stderr_pipe[0], buffer, BUFSIZE - 1)) > 0)
    {
      buffer[bytes_read] = '\0';
      callbackStdErr_(String(buffer));
    }

    close(stdout_pipe[0]);
    close(stderr_pipe[0]);

    // Status was already retrieved in the loop above

    if (WIFSIGNALED(status))
    {
      error_msg = "Process '" + exe + "' crashed hard (signal " + std::to_string(WTERMSIG(status)) + "). Please check the log.";
      if (verbose) callbackStdErr_(error_msg + '\n');
      return RETURNSTATE::CRASH;
    }
    else if (WIFEXITED(status))
    {
      int exit_code = WEXITSTATUS(status);
      if (exit_code == 127)
      {
        // Exit code 127 indicates command not found (exec failed)
        error_msg = "Process '" + exe + "' failed to start. Does it exist? Is it executable?";
        if (verbose) callbackStdErr_(error_msg + '\n');
        return RETURNSTATE::FAILED_TO_START;
      }
      else if (exit_code != 0)
      {
        error_msg = "Process '" + exe + "' did not finish successfully (exit code: " + std::to_string(exit_code) + "). Please check the log.";
        if (verbose) callbackStdErr_(error_msg + '\n');
        return RETURNSTATE::NONZERO_EXIT;
      }
    }
#endif

    if (verbose)
    {
      callbackStdOut_("Executed '" + exe + "' successfully!\n");
    }
    return RETURNSTATE::SUCCESS;
  }

} // namespace OpenMS
