import setuptools as st

with open("README.md", "r") as fh:
    long_description = fh.read()

st.setup(
    name="pyfocus",
    version="0.6",
    author="Nicholas Mancuso, Ruth Johnson",
    author_email="nicholas.mancuso@med.usc.com, ruthjohnson@ucla.com",
    description="Fine-map transcriptome-wide association studies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bogdanlab/focus",
    packages=st.find_packages(),
    package_data={'pyfocus': ['data/ld_blocks/*.bed']},
    install_requires=[
        # this is the minimum to perform fine-mapping given a prebuilt db.
        # functions that require addtl modules will warn/raise error
        # this is to minimize dependencies for the most-used scenario
        "opencv-python",
        "sqlalchemy",
        "matplotlib>=3.0.0",
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
