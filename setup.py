# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
exec(open("methylmap/version.py").read())

setup(
    name="methylmap",
    version=__version__,  # noqa: F821
    description="Plotting tool for population scale nucleotide modifications.",
    long_description=open(path.join(here, "README.md"), encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/EliseCoopman/methylmap",
    author="Elise Coopman",
    author_email="elisecoopman@yahoo.com",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    keywords="methylation plot",
    packages=find_packages() + ["scripts"],
    python_requires=">=3",
    install_requires=[
        "pandas",
        "numpy",
        "plotly>=5.4.0",
        "scipy",
        "dash",
        "pathlib",
        "tqdm",
        "dash_bootstrap_components",
    ],
    package_data={"methylmap": []},
    package_dir={"methylmap": "methylmap"},
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "methylmap=methylmap.methylmap:main",
        ],
    },
    data_files=[("", ["LICENSE"])],
)
