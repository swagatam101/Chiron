from setuptools import setup

setup(name='chiron',
      version='0.1.0',
      description='_CH_romatin _I_nformatics, _R_egional _O_rganization and _N_etwork',
      long_description = open('README.txt').read(),
      url='...',
      author='Swagatam Mukhopadhyay',
      author_email='swagatam.mukhopadhyay@gmail.com',
      license='MIT',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Toplogical analysis :: Chromatin Conformation Capture :: Analysis of Hi-C Data',
      ],
      keywords=['topology', 'chromatin', 'l1 sparsity', 'polymer', 'machine learning'], 
      packages=['chiron'],
      scripts=['bin/example_HiC.py'],
      install_requires=['numpy', 'matplotlib', 'networkx', 'scipy', 'pandas'],
      include_package_data=True,
      zip_safe=False)
      
      