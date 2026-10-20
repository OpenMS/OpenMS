// Registers the service worker of the demo application (sw.js, see there) with the directory of
// sw.js (examples/) as its scope. Call it before loading the module. Resolves to whether the worker
// controls the page; the page reloads once when it is not cross-origin isolated, so that the worker
// serves it with the isolation headers (threads).
export async function registerServiceWorker() {
  if (!("serviceWorker" in navigator) || location.protocol === "file:") return false;
  const RELOADED = "openms-service-worker-reloaded";
  try {
    // The scope is resolved against this module, as the worker's URL is. A relative scope alone
    // ("./") is resolved against the page, which gave each page directory a registration of its own.
    const scope = new URL("./", import.meta.url).href;
    await navigator.serviceWorker.register(new URL("sw.js", import.meta.url), { scope });
    await navigator.serviceWorker.ready;
    // Whether to reload depends on crossOriginIsolated, not on navigator.serviceWorker.controller:
    // the worker claims the open pages when it activates (clients.claim()), which can happen before
    // `ready` resolves. A claimed page has a controller, but the worker did not serve it, so it
    // lacks the isolation headers.
    if (crossOriginIsolated) {
      sessionStorage.removeItem(RELOADED);
    } else if (!sessionStorage.getItem(RELOADED)) {
      sessionStorage.setItem(RELOADED, "1");
      location.reload();
      await new Promise(() => {}); // the page is going away
    }
    return !!navigator.serviceWorker.controller;
  } catch (error) {
    console.warn(`OpenMS: the service worker could not be registered (${error.message}): the module is not cached, ` +
                 "and threads need the cross-origin isolation headers from the server.");
    return false;
  }
}
