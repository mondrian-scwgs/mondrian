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
            'breakpoint_utils = mondrian.utils.breakpoint_calling.utils:utils'
        ]
    },
    package_data={'': ['scripts/*.py', 'scripts/*.R', 'scripts/*.npz', "config/*.yaml", "data/*"]}
)
