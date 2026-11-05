# -*- coding: utf-8 -*-
"""
Run the usage examples in README.md and check the documented outputs, so that the README stays correct.

Every ```python block is executed in order in a shared namespace. For an expression statement followed by a comment
with a number or an ``array([...])``, the value is compared with the documented one at the documented precision.
"""
import ast
import os
import re

import numpy
import pytest

README = os.path.join(os.path.dirname(__file__), os.pardir, "README.md")


def _documented_values(comment):
    """Numbers in a documented output such as '3166.82...', 'array([ 934.83, 2969.91])' or '[nan, 2., 5.]'."""
    body = comment.split("    ")[0].replace("...", "")          # drop trailing remarks after the value
    if not re.match(r"^\s*(array\(|\[|-?\d|nan)", body):
        return None
    numbers = re.findall(r"nan|-?\d+\.?\d*(?:e-?\d+)?", body)
    if not numbers:
        return None
    decimals = max((len(n.split(".")[1]) if "." in n else 0) for n in numbers)
    return numpy.array([float(n) for n in numbers]), decimals


def _blocks():
    with open(README, encoding="utf-8") as f:
        return re.findall(r"```python\n(.*?)```", f.read(), re.S)


def test_readme_has_examples():
    assert len(_blocks()) >= 6


def test_readme_examples_run_and_match():
    namespace = {}
    checked = 0
    for block in _blocks():
        lines = block.splitlines()
        for node in ast.parse(block).body:
            if not isinstance(node, ast.Expr):
                exec(compile(ast.Module([node], []), "<README>", "exec"), namespace)
                continue
            value = eval(compile(ast.Expression(node.value), "<README>", "eval"), namespace)
            following = lines[node.end_lineno] if node.end_lineno < len(lines) else ""
            if not following.startswith("# "):
                continue
            documented = _documented_values(following[2:])
            if documented is None:
                continue
            expected, decimals = documented
            actual = numpy.atleast_1d(numpy.asarray(value, dtype=float))
            assert actual.shape == expected.shape, (ast.get_source_segment(block, node), actual, expected)
            numpy.testing.assert_allclose(actual, expected, atol=0.6 * 10.0 ** -decimals, equal_nan=True,
                                          err_msg=ast.get_source_segment(block, node))
            checked += 1
    assert checked >= 10
