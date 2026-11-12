// The open and the close of a host file in WasmFS' Node.js backend (mountHostDirectory()), linked with
// --js-library: Emscripten's own (src/lib/libwasmfs_node.js) return the errno of a failure as a positive
// number, while the backend (NodeState::open() and ::close() in system/lib/wasmfs/backends/node_backend.cpp)
// expects a negative one (https://github.com/emscripten-core/emscripten/issues/27825, Emscripten 6.0.12).
//
// - open: the backend took the errno of a failed open for a file descriptor of the Node.js process. A host
//   file that could not be opened, e.g. one deleted on the host after the module had looked it up (ENOENT,
//   44), was then open on the process' descriptor 44: reads went there, and closing the file closed that
//   descriptor, whatever it was. In the threaded module it was often one of the event loop of a thread's
//   Worker, and Node.js aborted at exit ("uv__io_poll: Assertion `errno == EINTR' failed").
// - close: a failing close of a host file (e.g. EBADF, 8) came out of FS.close() with a wrong errno (65528:
//   the backend negated the positive number, and WASI error codes are 16 bits).
//
// Host files of 2 GiB or more: the backend keeps the size of a host file and the offset and length of a read or write
// in 32 bits (uint32_t in node_backend.h, Emscripten 6.0.12; in wasm64 as well), and Emscripten's read and write took
// the offset for a signed number: a read at 2 GiB or beyond returned no bytes, as at the end of the file (a mounted
// mzML file of 3 GB was read as one of 2 GiB), and a write there failed. Here, read and write take offset and length as
// unsigned, which reaches 4 GiB, and the open refuses a host file of 4 GiB or more (EOVERFLOW, as a 32-bit program's
// open() without O_LARGEFILE), which the backend took for one of its size modulo 4 GiB, reading the wrong bytes. (A
// write at 4 GiB or beyond, into a host file that the module grows, still wraps around in the backend: offsets reach
// this file truncated to 32 bits.)
//
// Here, both return the negative errno, as the backend expects. This is stricter than Emscripten's, and
// intended:
// - FS.writeFile() of a read-only host file used to write to the stderr of the process (EACCES is 2 in
//   Emscripten, taken for the descriptor 2); it now traps (the abort in js_api.cpp: "Fatal error in
//   FS.writeFile").
// - open(O_CREAT|O_WRONLY, 0444) of a new host file used to succeed, on the descriptor 2 as well; it now
//   fails with EACCES.
//
// Drop this file (and its --js-library in CMakeLists.txt) once #27825 is fixed upstream. The check below
// reminds us to look whenever Emscripten is updated.

#if EMSCRIPTEN_VERSION != "6.0.12"
#warning wasmfs_node.js: Emscripten is not 6.0.12, recheck emscripten-core/emscripten#27825 (drop this file if fixed)
#endif

// A build without Node.js keeps Emscripten's stubs of these functions.
#if ENVIRONMENT_MAY_BE_NODE
addToLibrary({
  _wasmfs_node_open__deps: ["$wasmfsNodeConvertNodeCode"],
  _wasmfs_node_open: (path_p, mode_p) => {
    let fd;
    try {
      fd = fs.openSync(UTF8ToString(path_p), UTF8ToString(mode_p));
      if (fs.fstatSync(fd).size >= 2 ** 32) {
        fs.closeSync(fd);
        return -{{{ cDefs.EOVERFLOW }}};
      }
      return fd;
    } catch (e) {
      if (!e.code) throw e;
      if (fd !== undefined) {
        try { fs.closeSync(fd); } catch { /* the error of the open counts */ }
      }
      return -wasmfsNodeConvertNodeCode(e);
    }
  },

  _wasmfs_node_close__deps: ["$wasmfsNodeConvertNodeCode"],
  _wasmfs_node_close: (fd) => {
    try {
      fs.closeSync(fd);
    } catch (e) {
      if (!e.code) throw e;
      return -wasmfsNodeConvertNodeCode(e);
    }
    return 0;
  },

  _wasmfs_node_read__deps: ["$wasmfsTry"],
  _wasmfs_node_read: (fd, buf_p, len, pos, nread_p) => {
    return wasmfsTry(() => {
      let nread = fs.readSync(fd, HEAPU8, buf_p, len >>> 0, pos >>> 0);
      {{{ makeSetValue('nread_p', 0, 'nread', 'i32') }}};
      // implicitly return 0
    });
  },

  _wasmfs_node_write__deps: ["$wasmfsTry"],
  _wasmfs_node_write: (fd, buf_p, len, pos, nwritten_p) => {
    return wasmfsTry(() => {
      let nwritten = fs.writeSync(fd, HEAPU8, buf_p, len >>> 0, pos >>> 0);
      {{{ makeSetValue('nwritten_p', 0, 'nwritten', 'i32') }}};
      // implicitly return 0
    });
  },
});
#endif
