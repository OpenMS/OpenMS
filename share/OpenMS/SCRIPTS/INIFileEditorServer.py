#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# Local helper for the HTML INIFileEditor (share/OpenMS/HTML/INIFileEditor.html).
#
# Browsers without the File System Access API (Safari, Firefox) cannot write
# an opened file back to disk; their "Save" degrades to a download. This
# helper closes that gap for ALL browsers: it serves the editor page over
# http://127.0.0.1:<random port> and exposes the INI file given on the
# command line for reading and (in-place, atomic) writing. It also restores
# the command-line workflow of the Qt tool:
#
#     python3 INIFileEditorServer.py MyTool.ini        # opens default browser
#     python3 INIFileEditorServer.py --no-browser MyTool.ini   # print URL only
#
# Uses only the Python standard library. Stop with Ctrl+C.
#
# Security model:
#   - bound to 127.0.0.1; ephemeral port by default
#   - every /api request must carry a random per-session bearer token, which
#     is handed to the page in the URL fragment (never sent over the network,
#     never logged); this defeats CSRF from web pages and DNS rebinding
#   - Host and (when present) Origin headers are verified
#   - the server serves exactly one HTML file and exposes exactly the files
#     listed on the command line, addressed by index -- no paths ever cross
#     the API, so there is nothing to traverse
#   - writes are atomic (temp file + rename) and guarded by a modification
#     time check so concurrent external edits are detected, not clobbered
# --------------------------------------------------------------------------

import argparse
import hmac
import json
import os
import secrets
import sys
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path


def make_handler(html_bytes, files, token, allowed_hosts):
    class Handler(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"
        server_version = "INIFileEditorServer"
        sys_version = ""

        # ---- helpers -----------------------------------------------------
        MAX_BODY = 128 * 1024 * 1024

        def _deny(self, code, msg):
            body = msg.encode()
            self.send_response(code)
            self.send_header("Content-Type", "text/plain; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            # an unread request body would desync a kept-alive connection
            self.send_header("Connection", "close")
            self.close_connection = True
            self.end_headers()
            self.wfile.write(body)

        def _read_body(self):
            # drain the request body up front so error responses on kept-alive
            # connections never leave unread bytes in the stream
            if self.headers.get("Transfer-Encoding"):
                self._deny(411, "chunked request bodies are not supported; send Content-Length")
                return None
            try:
                length = int(self.headers.get("Content-Length", "0"))
            except ValueError:
                length = -1
            if length < 0 or length > self.MAX_BODY:
                self._deny(413, "request body missing a valid length or too large")
                return None
            return self.rfile.read(length)

        def _check_host_origin(self):
            host = self.headers.get("Host", "")
            if host not in allowed_hosts:
                self._deny(403, "forbidden (bad Host)")
                return False
            origin = self.headers.get("Origin")
            if origin is not None and origin not in ("http://" + h for h in allowed_hosts):
                self._deny(403, "forbidden (bad Origin)")
                return False
            return True

        def _check_token(self):
            got = self.headers.get("X-Auth-Token", "")
            if not hmac.compare_digest(got, token):
                self._deny(403, "forbidden (missing or wrong token)")
                return False
            return True

        def _send(self, code, body, ctype, extra=None):
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            for k, v in (extra or {}).items():
                self.send_header(k, v)
            self.end_headers()
            self.wfile.write(body)

        def _file_for_path(self):
            # /api/file/<index>
            parts = self.path.split("?", 1)[0].split("/")
            if len(parts) == 4 and parts[1] == "api" and parts[2] == "file" and parts[3].isdigit():
                idx = int(parts[3])
                if idx < len(files):
                    return idx, files[idx]
            return None, None

        # ---- HTTP methods ------------------------------------------------
        def do_GET(self):
            if not self._check_host_origin():
                return
            path = self.path.split("?", 1)[0]
            if path in ("/", "/INIFileEditor.html"):
                self._send(200, html_bytes, "text/html; charset=utf-8")
                return
            if path == "/api/files":
                if not self._check_token():
                    return
                listing = [{"id": i, "name": f.name, "mtime": str(f.stat().st_mtime_ns)}
                           for i, f in enumerate(files)]
                self._send(200, json.dumps({"files": listing}).encode(), "application/json")
                return
            idx, f = self._file_for_path()
            if f is not None:
                if not self._check_token():
                    return
                try:
                    data = f.read_bytes()
                except OSError as e:
                    self._deny(500, f"cannot read file: {e}")
                    return
                self._send(200, data, "application/xml",
                           {"X-File-Name": f.name, "X-File-Mtime": str(f.stat().st_mtime_ns)})
                return
            self._deny(404, "not found")

        def do_PUT(self):
            data = self._read_body()
            if data is None:
                return
            if not self._check_host_origin() or not self._check_token():
                return
            idx, f = self._file_for_path()
            if f is None:
                self._deny(404, "not found")
                return
            expected = self.headers.get("X-Expected-Mtime")
            try:
                current_mtime = str(f.stat().st_mtime_ns)
            except OSError as e:
                self._deny(500, f"cannot stat file: {e}")
                return
            if expected is not None and expected != current_mtime:
                self._send(409, json.dumps({"mtime": current_mtime}).encode(), "application/json")
                return
            try:
                tmp = f.with_name(f".{f.name}.tmp-{os.getpid()}")
                tmp.write_bytes(data)
                try:
                    os.chmod(tmp, f.stat().st_mode)
                except OSError:
                    pass
                with open(tmp, "rb") as fh:
                    os.fsync(fh.fileno())
                os.replace(tmp, f)
            except OSError as e:
                try:
                    tmp.unlink(missing_ok=True)
                except OSError:
                    pass
                self._deny(500, f"cannot write file: {e}")
                return
            new_mtime = str(f.stat().st_mtime_ns)
            self._send(200, json.dumps({"mtime": new_mtime}).encode(), "application/json",
                       {"X-File-Mtime": new_mtime})

        def do_POST(self):  # no POST endpoints
            self._deny(404, "not found")

        def log_message(self, fmt, *args):  # quiet; tokens must never hit a log
            pass

    return Handler


def main():
    parser = argparse.ArgumentParser(
        description="Serve the HTML INIFileEditor with in-place saving for any browser "
                    "(including Safari/Firefox, which lack the File System Access API).")
    parser.add_argument("ini", type=Path, help="INI file to edit")
    parser.add_argument("--port", type=int, default=0,
                        help="port to bind on 127.0.0.1 (default: random free port)")
    parser.add_argument("--html", type=Path, default=None,
                        help="path to INIFileEditor.html (default: ../HTML/INIFileEditor.html "
                             "relative to this script)")
    parser.add_argument("--no-browser", action="store_true",
                        help="do not open the default browser, only print the URL")
    args = parser.parse_args()

    ini = args.ini.resolve()
    if not ini.is_file():
        sys.exit(f"error: '{ini}' does not exist or is not a file")

    html_path = args.html or (Path(__file__).resolve().parent.parent / "HTML" / "INIFileEditor.html")
    if not html_path.is_file():
        sys.exit(f"error: editor page not found at '{html_path}' (use --html)")
    html_bytes = html_path.read_bytes()

    token = secrets.token_urlsafe(32)
    server = ThreadingHTTPServer(("127.0.0.1", args.port), None)
    port = server.server_address[1]
    allowed_hosts = {f"127.0.0.1:{port}", f"localhost:{port}"}
    server.RequestHandlerClass = make_handler(html_bytes, [ini], token, allowed_hosts)
    server.daemon_threads = True

    url = f"http://127.0.0.1:{port}/#t={token}"
    print(f"INIFileEditor helper serving '{ini}'")
    print(f"Editor URL: {url}")
    print("Open this URL in any browser (Safari/Firefox included); saving writes the file in place.")
    print("Press Ctrl+C to stop.")
    if not args.no_browser:
        threading.Timer(0.3, webbrowser.open, [url]).start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")


if __name__ == "__main__":
    main()
