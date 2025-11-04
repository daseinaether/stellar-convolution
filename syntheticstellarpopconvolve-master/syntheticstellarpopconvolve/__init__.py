import copy

from syntheticstellarpopconvolve.convolve import convolve  # noqa: F401
from syntheticstellarpopconvolve.default_convolution_config import (
    default_convolution_config as _default_convolution_config,
)
from syntheticstellarpopconvolve.default_convolution_instruction import (
    default_convolution_instruction as _default_convolution_instruction,
)

from ._version import __version__  # noqa: F401

default_convolution_config = copy.copy(_default_convolution_config)
default_convolution_instruction = copy.copy(_default_convolution_instruction)
