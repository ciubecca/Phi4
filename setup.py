from distutils.core import setup
from Cython.Build import cythonize

setup(
        ext_modules = cythonize("me.pyx", compiler_directives={"boundscheck":False})
        )
