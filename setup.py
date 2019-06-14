
from distutils.core import setup

setup(name='ad',
      description='Solve ODEs Using the Taylor Series Method (Automatic Differentiation)',
      author='Ian Smith',
      author_email='m4r35n357@gmail.com',
      version='1.0',
      url='https://github.com/m4r35n357/ODE-Playground',
      requires=['gmpy2', 'matplotlib', 'pillow', 'pi3d'],
      py_modules=['ad', 'playground'],
      )
