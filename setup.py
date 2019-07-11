import setuptools
import opencroplib

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='opencroplib',
    version=opencroplib.__version__,
    author='Jose A. Jimenez-Berni',
    author_email='berni@ias.csic.es',
    packages=['opencroplib'],
    url='https://github.com/OpenAgriTech/opencroplib',
    license='LICENSE.txt',
    description='Python Library for simulating physiological processes in plants',
    long_description=long_description,
    install_requires=[
        "python >= 2.7", 'numpy', 'pytz', 'pandas'
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
    ],
)
