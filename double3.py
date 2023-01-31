"""
Created on Fri Feb	4 20:05:59 2022

@author: Francesco Rizzuto
"""
# %%
import numpy as np

def scalar_prduct_for_array(a, b):
    """
    This function calculates the scalar product of two arrays using numpy's 
    sum function and passing the axis to sum along as 1.
    """
    return np.sum(a*b, axis=1)

def scalar_product(a, b):
    """
    This function calculates the scalar product of two objects of the Double3
    class by multiplying them together and returning the sum of their x, y, and z values.
    """
    c = a * b
    return c.x + c.y + c.z

def cross_product(a, b):
    """
    This function calculates the cross product of two objects of the Double3
    class by using the standard cross product formulas.
    """
    c = Double3()
    c.x = a.y * b.z - a.z * b.y
    c.y = a.z * b.x - a.x * b.z
    c.z = a.x * b.y - a.y * b.x
    return c

class Double3(object):
    def __init__(self, x=0.0, y=0.0, z=0.0):
        """
        This is the constructor for the Double3 class. It initializes the x, y, and z
        values to 0.0 by default, but they can also be set to specific values upon creation.
        """
        self.x = x
        self.y = y
        self.z = z
        
    def set_values(self, x=0.0, y=0.0, z=0.0):
        """
        This function allows the user to set the x, y, and z values of the Double3 object
        to specific values.
        """
        self.x = x
        self.y = y
        self.z = z 
        
    def as_array(self):   
        """
        This function returns the x, y, and z values of the Double3 object as a numpy array.
        """
        return np.array([self.x, self.y, self.z])
    
    #-------- SUM ---------    
    def __add__(self, other):
        """
        This function allows the user to add two Double3 objects together, or add a scalar
        value to a Double3 object.
        """
        result	 = Double3()
        if type(self) is type(other):
            result.x = self.x + other.x
            result.y = self.y + other.y
            result.z = self.z + other.z
        else:
            result.x = self.x + other
            result.y = self.y + other
            result.z = self.z + other
        return result

    def __radd__(self, other):
        """
        Overload the '+' operator for when the Double3 object is on the right side of the operator.
        """
        return self + other

    def __iadd__(self, other):
        """
        Overload the '+=' operator for the Double3 object.
        """
        if type(self) is type(other):
            self.x += other.x
            self.y += other.y
            self.z += other.z
        else:
            self.x += other
            self.y += other
            self.z += other
        return self

    def __sub__(self, other):
        """
        Overload the '-' operator for the Double3 object.
        """
        result = Double3()
        if type(self) is type(other):
            result.x = self.x - other.x
            result.y = self.y - other.y
            result.z = self.z - other.z
        else:
            result.x = self.x - other
            result.y = self.y - other
            result.z = self.z - other	 
        return result

    def __rsub__(self, other):
        """
        Overload the '-' operator for when the Double3 object is on the right side of the operator.
        """
        return self - other

    def __mul__(self, other):
        """
        Overload the '*' operator for the Double3 object.
        """
        result = Double3()
        if type(self) is type(other):
            result.x = self.x * other.x
            result.y = self.y * other.y
            result.z = self.z * other.z
        else:
            result.x = self.x * other
            result.y = self.y * other
            result.z = self.z * other
        return result

    def __rmul__(self, other):
        """
        Overload the '*' operator for when the Double3 object is on the right side of the operator.
        """
        return self * other

    def __truediv__(self, other):
        """
        Overload the '/' operator for the Double3 object.
        """
        result = Double3()
        result.x = self.x / other
        result.y = self.y / other
        result.z = self.z / other
        return result


    def __repr__(self):
        """
        Represent the Double3 object as a string in scientific notation.
        """
        return f'{self.x:.16e} {self.y:.16e} {self.z:.16e}'               


    def as_string(self):
        """return x, y, z vslues as single string

        Returns:
            string
        """

        return f'{self.x:.16e} {self.y:.16e} {self.z:.16e}'

# %%


