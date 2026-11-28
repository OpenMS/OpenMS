// Runs the workload in a dedicated Web Worker: the recommended way to use OpenMS in a page.
import { runWorkload } from "./workload.mjs";

onmessage = async ({ data }) => {
  try {
    postMessage({ result: await runWorkload(data) });
  } catch (e) {
    postMessage({ error: String(e?.stack ?? e) });
  }
};
