from setuptools  import setup, find_packages
import batracker

version_number = batracker.__version__

setup(name='batracker',
    version=version_number, 
    description='Open source Acoustic localisation of animal vocalisations',
    url='https://github.com/thejasvibr/batracker',
    author='Thejasvi Beleyur',
    license='MIT',
    packages=find_packages(),
    install_requires=['numpy', 'pandas','scipy','matplotlib','PyYAML'],
    zip_safe=False,
    include_package_data=True,
    classifiers=[ 'Intended Audience :: Science/Research',
        'Topic :: Multimedia :: Sound/Audio :: Analysis',
        'Topic :: Multimedia :: Sound/Audio',
		'Programming Language :: Python :: 3']
      )
