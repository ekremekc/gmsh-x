from setuptools import setup, find_packages

setup(
    name='gmsh_x',
    version = '0.0',
    author='Ekrem Ekici',
    author_email='ee331@cam.ac.uk',
    packages=['gmsh_x'],
    install_requires=[
        'gmsh',
        'h5py',
        'meshio',
        'numpy',
        'matplotlib',
        'scipy'
    ]
)
