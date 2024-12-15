from setuptools import setup

setup(
    name='mirimrs_destripe',
    version='0.1.0',    
    description='Small package to remove systematics based on non-zero background trace from JWST MIRI/MRS data',
    url='https://github.com/ecmatthews/mirimrs_destripe',
    author='Elisabeth Matthews',
    author_email='matthews@mpia.de',
    license='MIT license',
    packages=['mirimrs_destripe'],
    install_requires=['numpy',
                      'matplotlib',
                      'astropy',
                      'seaborn',
                      'jwst'
                      ],

)
