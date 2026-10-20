#!/usr/bin/env python3
"""Render a benchmark Markdown report and its local plots to PDF."""

import argparse
import html
import shutil
import subprocess
import tempfile
from pathlib import Path

import markdown


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = args.report.resolve()
    output = (args.output or report.with_suffix(".pdf")).resolve()
    browser = shutil.which("google-chrome") or shutil.which("chromium")
    if not browser:
        parser.error("PDF output needs Google Chrome or Chromium")
    body = markdown.markdown(report.read_text(encoding="utf-8"), extensions=["tables"])
    document = f"""<!doctype html><html lang="en"><meta charset="utf-8">
<title>OpenDIA benchmark report</title>
<base href="{html.escape(report.parent.as_uri() + '/')}">
<style>
@page {{ size: A4; margin: 14mm; }}
body {{ font: 10pt Arial, sans-serif; color: #222; }}
h1 {{ font-size: 20pt; }} h2 {{ font-size: 14pt; break-after: avoid; }}
p {{ line-height: 1.35; }} table {{ border-collapse: collapse; width: 100%; }}
th, td {{ border-bottom: 1px solid #ccc; padding: 5px; text-align: left; }}
img {{ width: 100%; max-height: 65mm; object-fit: contain; }}
table, img {{ break-inside: avoid; }} code {{ font-size: 8pt; }}
</style><body>{body}</body></html>"""
    # A private temporary profile avoids changes to the user's browser profile.
    with tempfile.TemporaryDirectory(prefix=".opendia-pdf-", dir=output.parent) as temporary:
        temporary = Path(temporary)
        page = temporary / "report.html"
        rendered = temporary / "report.pdf"
        page.write_text(document, encoding="utf-8")
        subprocess.run([
            browser, "--headless", "--disable-gpu", "--no-pdf-header-footer",
            f"--user-data-dir={temporary / 'profile'}",
            f"--print-to-pdf={rendered}", page.as_uri(),
        ], check=True, capture_output=True, text=True, timeout=60)
        if not rendered.is_file() or not rendered.read_bytes().startswith(b"%PDF-"):
            raise SystemExit("The browser did not create a PDF")
        rendered.replace(output)
    print(output)


if __name__ == "__main__":
    main()
