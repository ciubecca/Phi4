# File used to compile the Cython module

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

# extensions = [
    # Extension("me", ["me.pyx"], define_macros=[('CYTHON_TRACE', '1')])
# ]
# setup(
    # ext_modules = cythonize(extensions)
# )

setup(
    ext_modules = cythonize("me.pyx", compiler_directives={"boundscheck":False})
)




