#!/usr/bin/env python3
"""
Extract the full public API of pyopenms into a JSON file.

Usage:
    python extract_api.py output.json

This introspects every class, function, enum, and constant in pyopenms,
recording methods, properties, signatures, base classes, and enum values.
"""

import sys
import json
import inspect
import importlib
import enum


def get_signature(obj):
    """Try to extract a signature string for a callable."""
    try:
        sig = inspect.signature(obj)
        return str(sig)
    except (ValueError, TypeError):
        pass

    # For nanobind/pybind11/Cython objects, try __doc__ first line
    doc = getattr(obj, '__doc__', None)
    if doc:
        first_line = doc.strip().split('\n')[0]
        # Many C++ bindings put the signature in the first line of __doc__
        if '(' in first_line and ')' in first_line:
            return first_line
    return None


def get_overloaded_signatures(obj):
    """Extract all overloaded signatures from __doc__ (common in nanobind/pybind11)."""
    doc = getattr(obj, '__doc__', None)
    if not doc:
        return []

    sigs = []
    for line in doc.strip().split('\n'):
        line = line.strip()
        if '(' in line and ')' in line:
            # Check if it looks like a signature line (starts with function name or has ->)
            if '->' in line or line.endswith(')'):
                sigs.append(line)
            # Also catch lines that start with a number followed by a dot (nanobind overloads)
            elif line and line[0].isdigit() and '. ' in line:
                sigs.append(line.split('. ', 1)[1])
    return sigs


def extract_enum(cls):
    """Extract enum information."""
    info = {
        'type': 'enum',
        'values': {},
    }

    # Try standard Python enum
    if issubclass(type(cls), type) and issubclass(cls, enum.Enum):
        for member in cls:
            info['values'][member.name] = member.value
        return info

    # For C++ bound enums (pybind11/nanobind style), look for __members__
    members = getattr(cls, '__members__', None)
    if members:
        for name, val in members.items():
            try:
                info['values'][name] = int(val)
            except (TypeError, ValueError):
                info['values'][name] = str(val)
        return info

    # Try to get enum values from class attributes
    for name in sorted(dir(cls)):
        if name.startswith('_'):
            continue
        try:
            val = getattr(cls, name)
            if isinstance(val, (int, cls)):
                try:
                    info['values'][name] = int(val)
                except (TypeError, ValueError):
                    info['values'][name] = str(val)
        except Exception:
            pass

    return info


def is_enum_class(cls):
    """Check if a class is an enum (Python enum or C++ bound enum)."""
    if issubclass(type(cls), type) and issubclass(cls, enum.Enum):
        return True
    if hasattr(cls, '__members__'):
        return True
    # Heuristic: if all non-dunder attributes are instances of the class itself
    # (common pattern for pybind11/nanobind enums)
    attrs = [a for a in dir(cls) if not a.startswith('_')]
    if attrs:
        instances = sum(1 for a in attrs if isinstance(getattr(cls, a, None), cls))
        if instances > 0 and instances == len(attrs):
            return True
    return False


def extract_class(cls):
    """Extract class API information."""
    info = {
        'type': 'class',
        'bases': [],
        'methods': {},
        'properties': [],
        'static_methods': [],
        'class_methods': [],
        'nested_enums': {},
        'nested_classes': [],
    }

    # Base classes
    for base in getattr(cls, '__mro__', [cls])[1:]:
        base_name = base.__name__
        if base_name not in ('object', 'pybind11_object', 'instance'):
            info['bases'].append(f"{base.__module__}.{base_name}" if base.__module__ != 'builtins' else base_name)

    # Walk all attributes
    for name in sorted(dir(cls)):
        if name.startswith('__') and name.endswith('__'):
            # Only track special methods that matter for API compat
            if name in ('__init__', '__repr__', '__str__', '__len__', '__iter__',
                        '__getitem__', '__setitem__', '__contains__', '__eq__',
                        '__ne__', '__lt__', '__le__', '__gt__', '__ge__',
                        '__hash__', '__bool__', '__add__', '__iadd__',
                        '__mul__', '__imul__', '__copy__', '__deepcopy__'):
                pass
            else:
                continue
        elif name.startswith('_'):
            continue

        try:
            attr = getattr(cls, name)
        except Exception:
            continue

        # Check for nested enum
        if isinstance(attr, type) and is_enum_class(attr):
            info['nested_enums'][name] = extract_enum(attr)
            continue

        # Check for nested class
        if isinstance(attr, type):
            info['nested_classes'].append(name)
            continue

        # Property
        if isinstance(attr, property) or type(attr).__name__ in ('getset_descriptor', 'member_descriptor'):
            info['properties'].append(name)
            continue

        # nanobind properties show up differently
        attr_type_name = type(attr).__name__
        if attr_type_name in ('nb_method', 'nb_func'):
            sig = get_signature(attr)
            overloads = get_overloaded_signatures(attr)
            info['methods'][name] = {
                'signature': sig,
                'overloads': overloads if overloads else None,
            }
            continue

        # Callable (method, static method, etc.)
        if callable(attr):
            # Check if it's a static method
            is_static = isinstance(inspect.getattr_static(cls, name, None), staticmethod)
            is_classmethod = isinstance(inspect.getattr_static(cls, name, None), classmethod)

            sig = get_signature(attr)
            overloads = get_overloaded_signatures(attr)
            method_info = {
                'signature': sig,
                'overloads': overloads if overloads else None,
            }

            if is_static:
                info['static_methods'].append(name)
            elif is_classmethod:
                info['class_methods'].append(name)

            info['methods'][name] = method_info
            continue

        # If it's a nanobind/pybind11 property-like descriptor
        if hasattr(type(attr), '__get__'):
            info['properties'].append(name)
            continue

    # Remove empty fields to keep JSON clean
    info = {k: v for k, v in info.items() if v}
    info['type'] = 'class'

    return info


def extract_callable_singleton(obj):
    """Extract API from a CallableSingleton wrapper (used for DB singletons)."""
    info = {
        'type': 'class',
        'methods': {},
        'static_methods': ['getInstance'],
    }

    # Get the actual instance to introspect methods
    try:
        inst = obj.getInstance()
    except Exception:
        return info

    for name in sorted(dir(inst)):
        if name.startswith('_'):
            continue
        try:
            attr = getattr(inst, name)
        except Exception:
            continue
        if callable(attr):
            sig = get_signature(attr)
            overloads = get_overloaded_signatures(attr)
            info['methods'][name] = {
                'signature': sig,
                'overloads': overloads if overloads else None,
            }

    # getInstance itself
    info['methods']['getInstance'] = {'signature': '() -> instance', 'overloads': None}

    return info


def extract_api():
    """Extract the full pyopenms API."""
    import pyopenms

    api = {
        '_meta': {
            'version': getattr(pyopenms, '__version__', 'unknown'),
            'python': sys.version,
        },
        'classes': {},
        'functions': {},
        'constants': {},
        'enums': {},
    }

    for name in sorted(dir(pyopenms)):
        if name.startswith('_'):
            continue

        try:
            obj = getattr(pyopenms, name)
        except Exception as e:
            api['constants'][name] = f'<error: {e}>'
            continue

        # CallableSingleton wrapper (DB singletons like ElementDB, ResidueDB, etc.)
        if type(obj).__name__ == 'CallableSingleton':
            api['classes'][name] = extract_callable_singleton(obj)
            continue

        # Enum class
        if isinstance(obj, type) and is_enum_class(obj):
            api['enums'][name] = extract_enum(obj)
            continue

        # Regular class
        if isinstance(obj, type):
            api['classes'][name] = extract_class(obj)
            continue

        # Module (skip submodules)
        if inspect.ismodule(obj):
            continue

        # Function
        if callable(obj):
            sig = get_signature(obj)
            overloads = get_overloaded_signatures(obj)
            api['functions'][name] = {
                'signature': sig,
                'overloads': overloads if overloads else None,
            }
            continue

        # Constant
        api['constants'][name] = str(obj)[:200]

    return api


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} output.json", file=sys.stderr)
        sys.exit(1)

    output_path = sys.argv[1]

    print(f"Extracting pyopenms API...", file=sys.stderr)
    api = extract_api()

    # Summary
    print(f"  Version: {api['_meta']['version']}", file=sys.stderr)
    print(f"  Classes: {len(api['classes'])}", file=sys.stderr)
    print(f"  Enums: {len(api['enums'])}", file=sys.stderr)
    print(f"  Functions: {len(api['functions'])}", file=sys.stderr)
    print(f"  Constants: {len(api['constants'])}", file=sys.stderr)

    with open(output_path, 'w') as f:
        json.dump(api, f, indent=2, default=str)

    print(f"API written to {output_path}", file=sys.stderr)


if __name__ == '__main__':
    main()
