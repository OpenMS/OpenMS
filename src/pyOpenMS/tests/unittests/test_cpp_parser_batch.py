"""
Tests for libclang-based C++ header parsing performance features.

These focus on correctness of the batching/caching logic (not absolute timing).
"""

from __future__ import annotations

from pathlib import Path

import pytest

try:
    from generator.cpp_parser import CLANG_AVAILABLE
except Exception:
    CLANG_AVAILABLE = False


@pytest.mark.skipif(not CLANG_AVAILABLE, reason="libclang python bindings not available")
class TestCppHeaderParserBatching:
    def test_batch_parse_and_cache_roundtrip(self, tmp_path: Path):
        from generator.cpp_parser import CppHeaderParser

        header_dir = tmp_path / "headers"
        header_dir.mkdir(parents=True, exist_ok=True)

        a_h = header_dir / "A.h"
        b_h = header_dir / "B.h"

        a_h.write_text(
            """
namespace OpenMS
{
  class A
  {
  public:
    int foo() const;
  };
}
""".lstrip()
        )

        b_h.write_text(
            """
#include "A.h"

namespace OpenMS
{
  class B
  {
  public:
    void bar(const A& a);
  };
}
""".lstrip()
        )

        cache_dir = tmp_path / "libclang_cache"

        parser = CppHeaderParser(cache_dir=cache_dir)
        classes = parser.parse_headers(
            [a_h, b_h],
            batch_mode="on",
            batch_size=2,
            auto_batch_min_headers=1,
        )

        assert "A" in classes
        assert "B" in classes

        a_cls = classes["A"]
        b_cls = classes["B"]

        foo = a_cls.find_method("foo")
        assert foo is not None
        assert foo.return_type == "int"
        assert foo.is_const

        bar = b_cls.find_method("bar")
        assert bar is not None
        assert len(bar.parameters) == 1
        assert "A" in bar.parameters[0].type_str

        # Cache files should be written for both headers
        cache_files = list(cache_dir.glob("*.json"))
        assert len(cache_files) >= 2

        # Second run should be able to load purely from cache (no parsing)
        parser2 = CppHeaderParser(cache_dir=cache_dir)
        parser2._parse_headers_in_batch = lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("unexpected batch parse on cache hit")
        )
        parser2.parse_header = lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("unexpected per-header parse on cache hit")
        )

        classes2 = parser2.parse_headers(
            [a_h, b_h],
            batch_mode="on",
            batch_size=2,
            auto_batch_min_headers=1,
        )
        assert "A" in classes2
        assert "B" in classes2
