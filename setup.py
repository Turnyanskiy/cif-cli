from setuptools import setup
setup(
    name = 'cif_cli',
    version = '0.1.0',
    packages= ['cif_cli'],
    entry_points = {
        'console_scripts': [
            'pycli = pycli.__main__:main'
        ]
    }
)
