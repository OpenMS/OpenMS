// Service worker of the demo applications (examples/browser, examples/prose): each registers it
// with its own directory as scope (sw-client.mjs). It does two things a static web server
// cannot:
//
//  1. Threads. The threaded module needs a cross-origin isolated page, i.e. the response headers
//     Cross-Origin-Opener-Policy and Cross-Origin-Embedder-Policy on the page and on every
//     worker script, which hosts like GitHub Pages do not send. The worker adds them to those
//     responses (and to nothing else: see 2). On its first visit the page reloads once, so that
//     the worker serves it (sw-client.mjs).
//  2. Caching the module, and the models of a module with ONNX Runtime (models/*.onnx, fetched
//     when first used). openms.wasm is 30 MB, and Chrome keeps its compiled code only for a
//     module it could cache. Module files are served from the Cache API: on the first request
//     the worker stores the file and answers from the store at once, after which Chrome keeps the
//     compiled code and the module starts and runs faster from the second visit on (a search
//     against the human proteome: 3.7 s instead of 5.9 s). A response the worker alters (1) is not
//     eligible for that, so openms.wasm is passed on untouched; every other response gets the
//     headers (openms.mjs is also the script of the module's threads). The application's own
//     files are fetched from the network and fall back to the store offline.
//
// The stores are named after the package version (../package.json, read when the worker starts;
// an old version's stores are dropped). A hard reload (Shift+Reload) bypasses the worker.

const META_CACHE = "openms-meta";
const APP_DIR = new URL("./", self.location.href).pathname;            // examples/
const MODULE_FILE = /\/(openms\.(wasm|mjs|data)|index\.mjs|models\/[^/]+\.onnx)$/; // the module, the loader, its models

async function versionFromNetwork() {
  try {
    const response = await fetch(new URL("../package.json", self.location.href), { cache: "no-store" });
    return response.ok ? (await response.json()).version : null;
  } catch {
    return null; // offline
  }
}

/// The package version the stores are named after; drops the stores of an earlier one.
const version = (async () => {
  const meta = await caches.open(META_CACHE);
  const stored = await meta.match("version").then((r) => r?.text());
  const current = (await versionFromNetwork()) ?? stored ?? "0";
  if (current !== stored) {
    await meta.put("version", new Response(current));
    for (const name of await caches.keys()) {
      if (name !== META_CACHE && !name.endsWith(`-${current}`)) await caches.delete(name);
    }
  }
  return current;
})();

const moduleCache = async () => caches.open(`openms-module-${await version}`);
const appCache = async () => caches.open(`openms-app-${await version}`);

/// From the server, past the browser's HTTP cache, which may still hold the files of the previous
/// version (GitHub Pages: max-age=600): revalidated, i.e. a 304 where nothing changed.
const fromServer = (request) => fetch(request, { cache: "no-cache" });

/// From the store; the first time from the network, stored and answered from the store.
async function cacheFirst(request) {
  const cache = await moduleCache();
  const hit = await cache.match(request);
  if (hit) return hit;
  const response = await fromServer(request);
  if (!response.ok) return response;
  await cache.put(request, response);
  return (await cache.match(request)) ?? fromServer(request);
}

/// From the network, stored; from the store when the network fails.
async function networkFirst(request) {
  const cache = await appCache();
  try {
    const response = await fromServer(request);
    if (response.ok) await cache.put(request, response.clone());
    return response;
  } catch (error) {
    const hit = await cache.match(request);
    if (hit) return hit;
    throw error;
  }
}

function withIsolationHeaders(response) {
  if (response.status === 0) return response; // opaque
  const headers = new Headers(response.headers);
  headers.set("Cross-Origin-Embedder-Policy", "require-corp");
  headers.set("Cross-Origin-Opener-Policy", "same-origin");
  return new Response(response.body, { status: response.status, statusText: response.statusText, headers });
}

self.addEventListener("install", (event) => {
  event.waitUntil(version);
  self.skipWaiting();
});

self.addEventListener("activate", (event) => {
  event.waitUntil(self.clients.claim());
});

self.addEventListener("fetch", (event) => {
  const { request } = event;
  const url = new URL(request.url);
  if (request.method !== "GET" || url.origin !== self.location.origin) return;
  const isModule = MODULE_FILE.test(url.pathname);
  const isApp = url.pathname.startsWith(APP_DIR);
  if (!isModule && !isApp && request.mode !== "navigate") return;
  // (Chrome keeps the compiled code of openms.wasm only if its response is the stored one)
  const pristine = /\/openms\.(wasm|data)$/.test(url.pathname);
  const response = isModule ? cacheFirst(request) : isApp ? networkFirst(request) : fetch(request);
  event.respondWith(pristine ? response : response.then(withIsolationHeaders));
});
