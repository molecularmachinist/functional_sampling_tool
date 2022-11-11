from setuptools import setup, Extension
import numpy
from Cython.Build import cythonize


extension = Extension("functional_sampling_tool._ctransformations",
                      sources=["lib/_ctransformations.pyx"],
                      include_dirs=[numpy.get_include()],
                      extra_compile_args=["-O3"],
                      language="c++",
                      define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')])


setup(

    packages=["functional_sampling_tool",
              "functional_sampling_tool.analysis",
              "functional_sampling_tool.templates"],
    package_dir={"functional_sampling_tool": "src"},
    ext_modules=cythonize(extension, language_level="3")
)
