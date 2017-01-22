from distutils.core import setup
from Cython.Build import cythonize
# from Cython.Compiler.Options import directive_defaults

# directive_defaults['linetrace'] = True

setup(
    ext_modules = cythonize("oscillators.pyx")
)

