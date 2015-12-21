'''
HolyGrail (c) University of Manchester 2015

HolyGrail is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-few-public-methods
import sys

import holygrail.classification as classification
import synbiochem.optimisation.gen_alg as gen_alg


class ClassifierGeneticAlgorithm(gen_alg.GeneticAlgorithm):
    '''Class to optimise parameters for Classifier using a GA.'''

    def __init__(self, pop_size, sample_size, struct_patterns, split, args,
                 retain=0.2, random_select=0.05, mutate=0.01, verbose=False):
        '''Constructor.'''
        super(ClassifierGeneticAlgorithm, self).__init__(pop_size, args,
                                                         retain, random_select,
                                                         mutate, verbose)

        self.__classifier = classification.Classifier(sample_size,
                                                      struct_patterns,
                                                      split)

    def _fitness(self, individual):
        '''Determine the fitness of an individual.'''
        cls = self.__classifier.classify(**individual)

        if self._verbose:
            print str(cls[4]) + '\t' + str(individual)

        return 1 - cls[4]


def main(argv):
    '''main method.'''
    args = {'hidden_layers': range(0, 100),
            'learning_rate': [i * 0.001 for i in range(1, 10)],
            'momentum': [i * 0.1 for i in range(1, 10)],
            'patience': range(1, 10),
            'min_improvement': [i * 0.001 for i in range(1, 10)],
            'validate_every': range(1, 10),
            'batch_size': range(1, 10),
            'hidden_dropout': [i * 0.1 for i in range(0, 10)],
            'input_dropout': [i * 0.1 for i in range(0, 10)]}
    classifier = ClassifierGeneticAlgorithm(pop_size=int(argv[1]),
                                            sample_size=int(argv[2]),
                                            struct_patterns=argv[7:],
                                            split=float(argv[3]),
                                            args=args,
                                            retain=float(argv[4]),
                                            random_select=float(argv[5]),
                                            mutate=float(argv[6]),
                                            verbose=True)

    classifier.run()


if __name__ == '__main__':
    main(sys.argv)
