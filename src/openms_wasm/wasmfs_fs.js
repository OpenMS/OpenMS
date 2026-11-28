// FS.writeFile(), FS.write() and FS.read() of WasmFS (Emscripten's src/lib/libwasmfs.js), linked with --js-library.
//
// Emscripten's copy the data through a buffer from _malloc() and do not check it. The module grows its memory
// (ALLOW_MEMORY_GROWTH, hence ABORTING_MALLOC=0): where the memory cannot hold that buffer, e.g. a file larger than
// the memory has left, malloc returns 0, and the data was written from address 0 over the module's memory (FS.read():
// read from there). Here, they throw an FS.ErrnoError with ENOMEM (errno 48) instead, and the module stays usable.
// FS.write() also copies a typed array of bytes with one set instead of byte by byte (about 25 MB/s in Chrome), and
// FS.write() and FS.read() throw EINVAL (errno 28) for a range beyond their buffer.
//
// FS.write() and FS.read() are methods of the FS object, which this file cannot redefine: openmsWasmfsGuards replaces
// them when the module starts (__postset). FS_writeFile depends on it, so that it is linked whenever FS is.
// Drop this file (and its --js-library in CMakeLists.txt) once Emscripten checks these allocations. The check below
// reminds us to look whenever Emscripten is updated.

#if EMSCRIPTEN_VERSION != "6.0.12"
#warning wasmfs_fs.js: Emscripten is not 6.0.12, recheck the _malloc() of FS.writeFile/FS.write/FS.read in libwasmfs.js (drop this file if checked)
#endif

addToLibrary({
  $FS_writeFile__deps: ['$FS_fileDataToTypedArray', '_wasmfs_write_file', '$stackSave', '$stackRestore', '$stringToUTF8OnStack', 'malloc', 'free',
                        '$openmsWasmfsGuards'],
  $FS_writeFile: (path, data) => {
    data = FS_fileDataToTypedArray(data);
    var len = data.length;
    var dataBuffer = _malloc(len);
    if (!dataBuffer && len) {
      throw new FS.ErrnoError({{{ cDefs.ENOMEM }}});
    }
    var sp = stackSave();
    var pathBuffer = stringToUTF8OnStack(path);
    HEAPU8.set(data, dataBuffer);
    var ret = __wasmfs_write_file(pathBuffer, dataBuffer, len);
    _free(dataBuffer);
    stackRestore(sp);
    return ret;
  },

  $openmsWasmfsGuards__deps: ['$FS', 'malloc', 'free', '_wasmfs_read', '_wasmfs_pread', '_wasmfs_write', '_wasmfs_pwrite'],
  $openmsWasmfsGuards__postset: 'openmsWasmfsGuards();',
  $openmsWasmfsGuards: () => {
    // The range of the buffer is checked before the allocation, which a RangeError would otherwise leak (Node.js'
    // fs.readSync and fs.writeSync check it, too): EINVAL.
    var checkRange = (buffer, offset, length) => {
      if (!(offset >= 0 && length >= 0 && offset + length <= buffer.length)) {
        throw new FS.ErrnoError({{{ cDefs.EINVAL }}});
      }
    };
    FS.read = (stream, buffer, offset, length, position) => {
      checkRange(buffer, offset, length);
      var seeking = typeof position != 'undefined';
      var dataBuffer = _malloc(length);
      if (!dataBuffer && length) {
        throw new FS.ErrnoError({{{ cDefs.ENOMEM }}});
      }
      try {
        var bytesRead;
        if (seeking) {
          bytesRead = __wasmfs_pread(stream.fd, dataBuffer, length, {{{ splitI64('position') }}});
        } else {
          bytesRead = __wasmfs_read(stream.fd, dataBuffer, length);
        }
        if (bytesRead > 0) {
          buffer.set(HEAPU8.subarray(dataBuffer, dataBuffer + bytesRead), offset);
        }
      } finally {
        _free(dataBuffer);
      }
      return FS.handleError(bytesRead);
    };
    // canOwn is ignored, as by Emscripten's. A typed array of bytes is copied in one piece; other arrays (and typed
    // arrays of larger elements) value by value, as by Emscripten's.
    FS.write = (stream, buffer, offset, length, position, canOwn) => {
      checkRange(buffer, offset, length);
      var seeking = typeof position != 'undefined';
      var dataBuffer = _malloc(length);
      if (!dataBuffer && length) {
        throw new FS.ErrnoError({{{ cDefs.ENOMEM }}});
      }
      try {
        if (ArrayBuffer.isView(buffer) && buffer.BYTES_PER_ELEMENT === 1) {
          HEAPU8.set(new Uint8Array(buffer.buffer, buffer.byteOffset + offset, length), dataBuffer);
        } else {
          for (var i = 0; i < length; i++) {
            {{{ makeSetValue('dataBuffer', 'i', 'buffer[offset + i]', 'i8') }}};
          }
        }
        var bytesWritten;
        if (seeking) {
          bytesWritten = __wasmfs_pwrite(stream.fd, dataBuffer, length, {{{ splitI64('position') }}});
        } else {
          bytesWritten = __wasmfs_write(stream.fd, dataBuffer, length);
        }
      } finally {
        _free(dataBuffer);
      }
      return FS.handleError(bytesWritten);
    };
  },
});
