from distutils.core import setup, Extension

module = Extension ('factor',
                    libraries = ['gmp'],
                    sources=['factor.c'])

setup (name = 'Factor',
       version = '0.1',
       description = "A Python factorization function",
       author = "Eric Willisson",
       author_email = "ericwillisson@gmail.com",
       ext_modules = [module])
