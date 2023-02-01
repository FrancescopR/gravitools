from setuptools import find_packages, setup

setup(
    name='gravitools',
    version='1.0.0',
    author="Francesco Paolo Rizzuto",
    author_email="rizzutof93@gmail.com>",
    packages=find_packages(include=['gravitools/gravitools']),
    package_dir={'gravitools': 'gravitools/gravitools'},
    description='Python package to group a set of tools to study gravitational systems',
    license='MIT',
    install_requires=[],
#    setup_requires=['pytest-runner'],
#    tests_require=['pytest==4.4.1'],
    test_suite='tests',
)
