import versioneer
from setuptools import setup, find_packages

setup(
    name='mondrian',
    packages=find_packages(),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='single cell dna workflows',
    author='Diljot Grewal',
    author_email='diljot.grewal@gmail.com',
    entry_points={
        'console_scripts': [
            'mondrian = mondrian.run:main',
            'variant_utils = mondrian.utils.variant_calling.utils:utils',
            'breakpoint_utils = mondrian.utils.breakpoint_calling.utils:utils',
            'alignment_utils = mondrian.utils.alignment.utils:utils',
            'hmmcopy_utils = mondrian.utils.hmmcopy.utils:utils',
            'csverve_utils = mondrian.utils.io.csverve:utils',
            'pdf_utils = mondrian.utils.io.pdf:utils',
        ]
    },
    package_data={'': ['*.py', '*.R', '*.npz', "*.yaml", "data/*"]}
)
