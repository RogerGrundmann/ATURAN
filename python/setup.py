from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [
    Extension ("pyaturan",
              ['pyaturan.pyx', 'PythonStream.cpp'],
              language = 'c++',
              extra_compile_args = ["-std=c++11"],
              libraries = ['aturan'],
              include_dirs = ['../planet', '../lib', '../tinyxml2'],
              library_dirs = ['..'],
              extra_link_args = ['-fopenmp'],
              )]

setup(
    name = 'pyaturan',
    ext_modules = cythonize(extensions, language_level = "2"),
    )
