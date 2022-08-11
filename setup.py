#taken heavily from:
#   https://github.com/pybind/python_example
#   https://pybind11.readthedocs.io/en/stable/compiling.html#building-with-setuptools
from pathlib import Path
import sys

from pybind11 import get_cmake_dir
# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

__version__ = "0.0.1"
NUM_MODES = 8000
NUM_BITS = 8
include_dirs = Path('./ext/boost/')
# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)

ext_modules = [
    Pybind11Extension(
        name="pyket",
        sources=sorted(Path('.').glob('**/*.cpp')),
        # Example: passing in the version to the compiled code
        define_macros = [
            ('VERSION_INFO', __version__),
            ('NUM_MODES', NUM_MODES),
            ('NUM_BITS',NUM_BITS),
            ('MODULE_NAME',__name__)
            ],
        ),
        include_dirs = [
            Path('./extern/rapidjson/include'),
            Path('./extern/boost'),
            Path('./extern/eigen'),
        ]
]

setup(ext_modules=ext_modules,cmdclass={"build_ext": build_ext})
