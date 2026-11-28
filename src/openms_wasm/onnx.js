// A module with ONNX Runtime (WITH_ONNX), linked with --pre-js before threads.js.

/**
 * Threads of ONNX Runtime, started when OpenMS first uses it (e.g. ProSE's peptdeep:enable): one
 * pool for all sessions with as many threads as OpenMP computes with, the calling one included
 * (src/openms/source/ML/ONNX/ONNXEnvironment.cpp). threads.js adds their Web Workers to the pool.
 */
var openmsOnnxThreads = true;
