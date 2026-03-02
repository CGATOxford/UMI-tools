from setuptools import build_meta as _orig
from setuptools.build_meta import *


def get_cython_requirement(config_settings=None):
    """
    Parse a mapping containing configuration passed to the build system.
    This mapping can be populated by pip's '--config-settings' argument.
    e.g.: `--config-settings with-cython=true`
    """
    extra_requirements = []
    if config_settings and config_settings.get("with-cython") == "true":
        extra_requirements = ["cython"]
    return extra_requirements

def get_requires_for_build_wheel(config_settings=None):
    """
    The requirements for building a wheel, optionally with cython installed.
    Adding cython as build requirements will cause the .pyx source files 
    to be compiled to C code which is added to the wheel.
    """
    return _orig.get_requires_for_build_wheel(config_settings) \
        + get_cython_requirement(config_settings)


def get_requires_for_build_sdist(config_settings=None):
    """
    The requirements for building a source distribution (.tar.gz), 
    optionally with cython installed. Adding cython as build 
    requirements will cause the .pyx source files to be compiled 
    to C code which is added to the tar file.
    """
    return _orig.get_requires_for_build_sdist(config_settings) \
        + get_cython_requirement(config_settings)