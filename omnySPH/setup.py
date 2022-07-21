import sys

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

from setuptools import find_packages

setup(
    name="omnySPH",
    version="0.0.2",
    description="a minimal example package (with pybind11)",
    author="Henry Schreiner",
    license="MIT",
    packages=['omnySPH'],
    package_dir={"": "src"},
    cmake_install_dir="src/omnySPH",
    python_requires=">=3.6",
)
