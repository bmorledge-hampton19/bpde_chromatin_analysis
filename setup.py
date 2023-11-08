from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="bpde_chromatin_analysis",
    version="0.1",
    description='Analyze BPDE lesions with respect to genomic features like nucleosomes',
    long_description_content_type="text/markdown",
    url='https://github.com/bmorledge-hampton19/bpde_chromatin_analysis',
    author='Ben Morledge-Hampton',
    author_email='b.morledge-hampton@wsu.edu',
    license='MIT',
    python_requires='>=3.10',
    packages=find_packages(exclude=["bpde_chromatin_analysis._notebooks"]),
    install_requires=["benbiohelpers", "mutperiodpy"]
)