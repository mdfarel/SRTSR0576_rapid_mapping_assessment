"""Local preview of the web GIS (webgis/) with HTTP range requests, which PMTiles needs.

Python's built-in file server always sends whole files; PMTiles reads each tile archive in small
byte ranges, so this server answers a "Range: bytes=start-end" request with just those bytes.

    python serve.py [port]        (or: make serve), then open http://localhost:8765
"""
import os
import re
import sys
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

WEB = Path(__file__).resolve().parent / "webgis"


class RangeHandler(SimpleHTTPRequestHandler):
    def send_head(self):
        """Headers for a 206 partial response when the request names a byte range; otherwise the usual reply."""
        path = self.translate_path(self.path)
        m = re.match(r"bytes=(\d+)-(\d*)", self.headers.get("Range", ""))
        if not m or not os.path.isfile(path):
            return super().send_head()
        size = os.path.getsize(path)
        start = int(m.group(1))
        end = min(int(m.group(2)) if m.group(2) else size - 1, size - 1)
        if start >= size:
            self.send_error(416)
            return None
        f = open(path, "rb")
        f.seek(start)
        self.send_response(206)
        self.send_header("Content-Type", self.guess_type(path))
        self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.send_header("Content-Length", str(end - start + 1))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()
        self._remaining = end - start + 1
        return f

    def copyfile(self, source, outputfile):
        """Send only the requested bytes after a range request, and the whole file otherwise."""
        remaining = getattr(self, "_remaining", None)
        if remaining is None:
            return super().copyfile(source, outputfile)
        while remaining > 0:
            chunk = source.read(min(65536, remaining))
            if not chunk:
                break
            outputfile.write(chunk)
            remaining -= len(chunk)
        self._remaining = None


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8765
    print(f"http://localhost:{port}")
    ThreadingHTTPServer(("", port), partial(RangeHandler, directory=str(WEB))).serve_forever()
