from setuptools import setup, find_packages

setup(
    name='ncbi',
    version='0.1.1',
    packages= find_packages(),
    include_package_data=True,
    install_requires=[
        'Click', "biopython", "requests",
        "pandas", "beautifulsoup4", "plotly",
        "fpdf", "kaleido"
    ],
    entry_points={
        'console_scripts': [
            'ncbi = ncbi.main:cli',
        ],
    },
)
