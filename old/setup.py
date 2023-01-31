from setuptools import find_packages, setup
setup(
    name='gravitools',
    packages=find_packages(include=['gravitools']),
    version='0.1.0',
    description='Python package to performe and analyze gravitational scattering experiments',
    author='Francesco Paolo Rizzuto',
    license='MIT',
    install_requires=[],
    setup_requires=['pytest-runner'],
    tests_require=['pytest==4.4.1'],
    test_suite='tests',
)
