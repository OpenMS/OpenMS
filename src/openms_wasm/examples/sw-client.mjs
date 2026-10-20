// Registers the service worker of the demo application (sw.js, see there) for the calling
// page's directory. Call it before loading the module. Resolves to whether the worker serves
// the page; the page reloads once when it needs the worker's isolation headers (threads).
export async function registerServiceWorker() {
  if (!("serviceWorker" in navigator) || location.protocol === "file:") return false;
  const RELOADED = "openms-service-worker-reloaded";
  try {
    await navigator.serviceWorker.register(new URL("sw.js", import.meta.url), { scope: "./" });
    await navigator.serviceWorker.ready;
    if (navigator.serviceWorker.controller) {
      sessionStorage.removeItem(RELOADED);
      return true;
    }
    if (!crossOriginIsolated && !sessionStorage.getItem(RELOADED)) {
      sessionStorage.setItem(RELOADED, "1");
      location.reload();
      await new Promise(() => {}); // the page is going away
    }
    return false;
  } catch (error) {
    console.warn(`OpenMS: the service worker could not be registered (${error.message}): the module is not cached, ` +
                 "and threads need the cross-origin isolation headers from the server.");
    return false;
  }
}
