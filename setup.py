from setuptools import setup

def get_version_number():
    main_ns = {}
    for line in open('pyrasorbitalkit/__init__.py', 'r'):
        if (line.startswith('__version__')):
            exec(line, main_ns)
            return main_ns['__version__']

#Make python package
setup(
    name='pyrasorbitalkit',
    version=get_version_number(),
    description="Python package for different RASCI calculations with NOs",
    long_description=open('README.md').read(),
    install_requires=['numpy', 'scipy', 'pyqchem'],
    author='Aaron Rodriguez',
    author_email='aarodjim@gmail.com',
    packages=['pyrasorbitalkit',
              ]
)

