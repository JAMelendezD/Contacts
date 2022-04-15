import unittest
import numpy as np
from package import distance_mat, contacts

class TestDistanceMat(unittest.TestCase):
    def test_inputs(self):
        self.assertRaises(TypeError, distance_mat, [1,0,0], 2)
        self.assertRaises(TypeError, distance_mat, 5, 5)
        self.assertRaises(TypeError, distance_mat, 3+5j, 5)
        self.assertRaises(TypeError, distance_mat, np.array([[1,0,0]]), -1)
        self.assertRaises(TypeError, distance_mat, np.array([[1,0,0]]), 3+5j)
    def test_valid_input(self):
        self.assertRaises(ValueError, distance_mat, np.array([[1,0,0]]), 2)
    def test_answer(self):
        np.testing.assert_allclose(
            distance_mat(np.array([[0,0,0],[0,1,0],[0,0,1]]),5),
            np.array([[0.,1.,1.],[1.,0.,1.41421356],[1.,1.41421356,0.]]))
        np.testing.assert_allclose(
            distance_mat(np.array([[0,0],[0,9]]),10),
            np.array([[0.,1.],[1.,0.]]))
        np.testing.assert_allclose(
            distance_mat(np.array([[0,0,0],[0,0,4]]),5),
            np.array([[0.,1.],[1.,0.]]))


class TestContacts(unittest.TestCase):
    def test_answer(self):
        np.testing.assert_allclose(
            contacts(np.array([[0,0,0]]),np.array([[0,0,1]]),1),
            np.array([[1]]))
        np.testing.assert_allclose(
            contacts(np.array([[0,0,0]]),np.array([[0,0,1]]),0.9),
            np.array([[0]]))
        np.testing.assert_allclose(
            contacts(np.array([[0,0,0],[1,0,0],[0,2,0]]),
            np.array([[0,0,1],[1,0,0],[0,0,0]]),1),
            np.array([[1,1,1],[0,1,1],[0,0,0]]))
        np.testing.assert_allclose(
            contacts(np.array([[0,0,0],[1,0,0],[0,2,0]]),
            np.array([[0,0,1],[1,0,0],[0,0,0]]),5),
            np.array([[1,1,1],[1,1,1],[1,1,1]]))
        np.testing.assert_allclose(
            contacts(np.array([[0,0,0],[1,0,0],[0,2,0]]),
            np.array([[0,0,1],[1,0,0],[0,0,0]]),0.0),
            np.array([[0,0,1],[0,1,0],[0,0,0]]))