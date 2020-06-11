
"""A module for creating new functions from a strings.

Index
-----
.. currentmodule:: CAT.data_handling.str_to_func
.. autosummary::
    str_to_func
    compile_func

API
---
.. autofunction:: str_to_func
.. autofunction:: compile_func

"""
from types import FunctionType, CodeType
from typing import Union, Callable

import numpy as np

from .warn_map import WARN_MAP

__all__ = ['str_to_func']


def str_to_func(func_str: Union[Callable[[np.ndarray], np.ndarray], str],
                shape: str = 'raise', dtype: str = 'raise') -> Callable[[np.ndarray], np.ndarray]:
    """Convert a string-representation of a function into an actual function and validate its output."""  # noqa
    shape_func = WARN_MAP[shape]
    dtype_func = WARN_MAP[dtype]

    # Create the new function
    if callable(func_str):  # **func_str** is already a function; no further compilation required
        func = func_str
    else:
        func = compile_func(func_str)

    # Validate and return the function
    ar1 = np.array([[0.0, np.nan], [np.inf, 5.0]])
    ar2 = func(ar1)
    if ar1.shape != ar2.shape:
        shape_func(ValueError(f"{func_str} -> y: input and output arrays have "
                              f"non-identical shapes: {ar1.shape} and {ar2.shape}"))
    if ar1.dtype != ar2.dtype:
        dtype_func(ValueError(f"{func_str} -> y: input and output arrays have non-identical "
                              f"data types: '{ar1.dtype.name}' and '{ar2.dtype.name}'"))
    return func


def compile_func(func_str: str) -> FunctionType:
    """Compile **func_str** into a :class:`code<types.CodeType>` instance."""
    global_dict = {'np': np}
    x = np.array([[0.0, np.nan], [np.inf, 5.0]])

    func_str = func_str.replace('numpy', 'np')

    # Load all module specified in **func_str**
    while True:
        try:
            x.copy()
            exec(func_str)
        except NameError as ex:
            module = str(ex).split("'")[1]
            exec(f'import {module}')
            global_dict[module] = eval(module)
        else:
            break

    # Exchange all ';' characters for '\n'
    str_list = func_str.split(';')
    func_str2 = ''.join(f'    {i.strip().rstrip()}\n' for i in str_list[:-1])
    func_str2 += f'    return {str_list[-1].strip().rstrip()}'

    # Compile the code of the to-be returned function
    code_compile = compile(f'def weight(x):\n{func_str2}', '<string>', 'exec')
    for code in code_compile.co_consts:
        if isinstance(code, CodeType):
            break
    else:
        raise ValueError("Failed to find 'code' instance in 'code_compile.co_consts'")

    # Construct and return the new function
    func = FunctionType(code, globals=global_dict, name='weight')
    func.__doc__ = f"Dynamically generated function; returns ``{func_str}``."
    func.__annotations__ = {'x': np.ndarray, 'return': np.ndarray}
    return func
