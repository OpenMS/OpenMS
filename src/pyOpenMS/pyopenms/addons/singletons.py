"""Singleton wrapper for DB types that use getInstance()."""

from __future__ import annotations


def make_singleton_callable(original_cls):
    """Wrap a singleton class so cls() returns cls.getInstance().

    Uses a callable wrapper that delegates attribute access to the original class,
    and calling it returns the singleton instance.
    """
    class CallableSingletonMeta(type):
        """Metaclass that makes isinstance() work against the original class."""

        def __instancecheck__(cls, instance):
            return isinstance(instance, original_cls)

    class CallableSingleton(metaclass=CallableSingletonMeta):
        """Makes ClassName() return ClassName.getInstance()."""

        def __call__(self):
            return original_cls.getInstance()

        def __getattr__(self, name):
            return getattr(original_cls, name)

        def __repr__(self):
            return repr(original_cls)

    wrapper = CallableSingleton()
    # lets introspection tools (the API documentation) reach the wrapped class; not
    # __wrapped__, which inspect.signature() would follow to the nanobind class
    wrapper._singleton_class = original_cls
    return wrapper


SINGLETON_CLASSES = [
    "ElementDB", "RibonucleotideDB", "ModificationsDB",
    "ProteaseDB", "RNaseDB", "CrossLinksDB",
    "EnzymesDB", "ResidueDB",
]
