import unittest
import numpy as np
from package import distance

class TestDistance(unittest.TestCase):
    def test_inputs(self):
        self.assertRaises(TypeError, distance, [1,0,0], 2)
        self.assertRaises(TypeError, distance, 5, 5)
        self.assertRaises(TypeError, distance, 3+5j, 5)
        self.assertRaises(TypeError, distance, np.array([[1,0,0]]), -1)
        self.assertRaises(TypeError, distance, np.array([[1,0,0]]), 3+5j)
    def test_valid_input(self):
        self.assertRaises(ValueError, distance, np.array([[1,0,0]]), 2)
    def test_answer(self):
        np.testing.assert_allclose(
            distance(np.array([[0,0,0],[0,1,0],[0,0,1]]),5),
            np.array([[0.,1.,1.],[1.,0.,1.41421356],[1.,1.41421356,0.]]))
        np.testing.assert_allclose(
            distance(np.array([[0,0],[0,9]]),10),
            np.array([[0.,1.],[1.,0.]]))
        np.testing.assert_allclose(
            distance(np.array([[0,0,0],[0,0,4]]),5),
            np.array([[0.,1.],[1.,0.]]))
