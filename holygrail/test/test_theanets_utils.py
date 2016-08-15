'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-many-public-methods
import random
import unittest

from sklearn import datasets
from sklearn.datasets.samples_generator import make_blobs
import numpy

import holygrail.theanets_utils as theanets_utils


class Test(unittest.TestCase):
    '''Tests the ann module.'''

    def test_randomise_order(self):
        '''Tests the randomise_order method.'''
        list_a = [int(random.random() * 1000) for _ in range(0, 1024)]
        list_b = [v for v in list_a]
        list_a, list_b = theanets_utils.randomise_order(list_a, list_b)

        for i, val in enumerate(list_a):
            self.assertEqual(val, list_b[i])


class TestClassifier(unittest.TestCase):
    '''Tests the Classifier class.'''

    def test_classify(self):
        '''Tests the classify method.'''

        # Generate some data, convert to list of floats (inputs) and string
        # 'names' for output classifications:
        x_data, y_data = make_blobs(n_samples=1000, centers=5, n_features=3,
                                    cluster_std=1.0, random_state=0)
        x_data = x_data.tolist()
        y_data = [str(unichr(y + ord('A'))) for y in y_data]

        # Split data into training and classifying:
        ind = int(0.8 * len(x_data))

        classifier = theanets_utils.Classifier(x_data[:ind], y_data[:ind])
        classifier.train(hyperparams={'optimize': 'sgd',
                                      'learning_rate': 5e-4})

        y_test = y_data[ind:]
        y_pred, _, _, _, _ = classifier.classify(x_data[ind:], y_test)
        self.assertTrue(sum([i == j for i, j in zip(y_test, y_pred)]) /
                        float(len(y_test)) > 0.9)


class TestRegressor(unittest.TestCase):
    '''Tests the Regressor class.'''

    def test_predict(self):
        '''Tests the predict method.'''

        # Load the diabetes dataset:
        dataset = datasets.load_diabetes()

        x_data = dataset.data.tolist()
        y_data = dataset.target.tolist()

        x_data, y_data = theanets_utils.randomise_order(x_data, y_data)

        # Split data into training and classifying:
        ind = int(0.8 * len(x_data))

        y_train = [[y] for y in y_data[:ind]]
        regressor = theanets_utils.Regressor(x_data[:ind], y_train)

        regressor.train(hidden_layers=[1])
        y_pred = regressor.predict(x_data[ind:])

        self.assertTrue(numpy.sqrt(numpy.mean((y_data[ind:] - y_pred) ** 2)) <
                        200)
