// The open of a host file in WasmFS' Node.js backend (mountHostDirectory()), linked with --js-library:
// Emscripten's own (src/lib/libwasmfs_node.js) returns the errno of a failed open as a positive
// number, which the backend (NodeState::open() in system/lib/wasmfs/backends/node_backend.cpp) takes
// for a file descriptor (https://github.com/emscripten-core/emscripten/issues/27825, Emscripten 6.0.12).
// A host file that could not be opened, e.g. one deleted on the host after the module had looked it up
// (ENOENT, 44), was then open on the Node.js process' descriptor 44: reads went there, and closing the
// file closed that descriptor, whatever it was. In the threaded module it was often one of the event
// loop of a thread's Worker, and Node.js aborted at exit ("uv__io_poll: Assertion `errno == EINTR'
// failed"). Here, a failed open returns the negative errno, as the backend expects.

addToLibrary({
  _wasmfs_node_open__deps: ["$wasmfsNodeConvertNodeCode"],
  _wasmfs_node_open: (path_p, mode_p) => {
    try {
      return fs.openSync(UTF8ToString(path_p), UTF8ToString(mode_p));
    } catch (e) {
      if (!e.code) throw e;
      return -wasmfsNodeConvertNodeCode(e);
    }
  },
});
