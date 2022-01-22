import setuptools as st

with open("README.md", "r") as fh:
    long_description = fh.read()


st.setup(
    name="pyfocus",
    version="0.8",
    author="Nicholas Mancuso, Ruth Johnson, Zeyun Lu",
    author_email="nicholas.mancuso@med.usc.com, ruthjohnson@ucla.com, zeyunlu@usc.edu",
    description="Fine-map transcriptome-wide association studies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mancusolab/focus",
    packages=st.find_packages(),
    package_data={'pyfocus': ['data/ld_blocks/*.bed', "data/gencode_map_v37.tsv"]},
    install_requires=[
        # this is the minimum to perform fine-mapping given a prebuilt db.
        # functions that require addtl modules will warn/raise error
        # this is to minimize dependencies for the most-used scenario
        "opencv-python",
        "sqlalchemy",
        "matplotlib>=3.1.0",
        "seaborn",
        "numpy",
        "scipy",
        "pandas>=0.23.0",
        "pandas-plink"
      ],
    scripts=[
        "bin/focus",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
