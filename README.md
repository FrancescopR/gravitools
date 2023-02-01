Editable install from local folder

pip install -e .




This Package is a simple example to illustrate how to create a python package.
(see https://medium.com/analytics-vidhya/how-to-create-a-python-library-7d5aea80cc3f )

1) prerequisites

pip install wheel
pip install setuptools
pip install twine


2) create a root folder (wat is the name of this folder? )

3) content

- Create an empty file called setup.py
- reate an empty file called README.md
- Create a folder called "nameofyourpackage"
- Create an empty file inside "nameofyourpackage" that is called __init__.py
- Also, in the same folder, create all your modueles and files .py
- create a folder tests in your root folder. Inside, create an empty __init__.py file and an empty test\_myfunctions.py


4) installing packages for testing

- pip install pytest==4.4.1
- pip install pytest-runner==4.4

5) Run test 
- python setup.py pytest

6) build the library
- python setup.py bdist\_wheel

7) install the library
- pip install dist/addnumbers-0.1.0-py3-none-any.whl
