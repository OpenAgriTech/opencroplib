from distutils.core import setup
import opencroplib

setup(
    name='opencroplib',
    version=opencroplib.__version__,
    author='Jose A. Jimenez-Berni',
    author_email='berni@ias.csic.es',
    packages=['opencroplib'],
    url='http://pypi.python.org/pypi/opencroplib/',
    license='LICENSE.txt',
    description='Python Library for simulating physiological processes in plants',
    long_description=open('README.md').read(),
    install_requires=[
        "python >= 2.7", 'numpy', 'pytz', 'pandas'
    ],
)
