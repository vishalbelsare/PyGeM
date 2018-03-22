from setuptools import setup, find_packages


def readme():
    """
    This function just return the content of README.md
    """
    with open('README.md') as f:
        return f.read()


setup(
    name='pygem',
    version='0.2',
    description='Tools for gemetrical morphing.',
    long_description=readme(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics'
    ],
    keywords='dimension_reduction mathematics ffd morphing iges stl vtk openfoam',
    url='https://github.com/mathLab/PyGeM',
    author='Marco Tezzele, Nicola Demo',
    author_email='marcotez@gmail.com, demo.nicola@gmail.com',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'numpy', 'numpy-stl', 'scipy', 'matplotlib', 'vtk', 'enum34',
        'Sphinx==1.4', 'sphinx_rtd_theme'
    ],
    test_suite='nose.collector',
    tests_require=['nose'],
    include_package_data=True,
    zip_safe=False)
