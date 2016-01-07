'''
HolyGrail (c) University of Manchester 2015

HolyGrail is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-few-public-methods
import sys

import holygrail.classification as classification
import holygrail.data
import synbiochem.optimisation.gen_alg as gen_alg


class ClassifierGeneticAlgorithm(gen_alg.GeneticAlgorithm):
    '''Class to optimise parameters for Classifier using a GA.'''

    def __init__(self, pop_size, pdb_data, split, args,
                 retain=0.2, random_select=0.05, mutate=0.01, verbose=False):
        '''Constructor.'''
        super(ClassifierGeneticAlgorithm, self).__init__(pop_size, args,
                                                         retain, random_select,
                                                         mutate, verbose)

        self.__classifier = classification.Classifier(pdb_data, split)

    def _fitness(self, individual):
        '''Determine the fitness of an individual.'''

        # Form hidden layers array:
        num_hidden_layers = individual.pop('num_hidden_layers', 1)
        activ_func = individual.pop('activ_func', 'relu')
        num_nodes = individual.pop('num_nodes', 100)
        individual['hidden_layers'] = [(num_nodes, activ_func)] * \
            num_hidden_layers

        cls = self.__classifier.classify(**individual)

        # Reform individual dict:
        individual.pop('hidden_layers')
        individual['num_hidden_layers'] = num_hidden_layers
        individual['activ_func'] = activ_func
        individual['num_nodes'] = num_nodes

        if self._verbose:
            print str(cls[4]) + '\t' + str(individual)

        return 1 - cls[4]


def main(argv):
    '''main method.'''

    # Get random peptides that match structure patterns from PDB:
    pdb_data, hammings = holygrail.data.sample_seqs(int(argv[2]),
                                                    argv[7:])

    print hammings

    args = {'aa_props_filter': range(1, (2**holygrail.NUM_AA_PROPS)),
            'input_noise': [i * 0.1 for i in range(0, 10)],
            'hidden_noise': [i * 0.1 for i in range(0, 10)],
            'num_hidden_layers': range(1, 4),
            'num_nodes': range(1, 100),
            'activ_func': ['linear', 'sigmoid', 'logistic', 'tanh', 'softplus',
                           'softmax', 'relu', 'rect:min', 'rect:max',
                           'norm:mean', 'norm:max', 'norm:std', 'norm:z',
                           'prelu', 'lgrelu'],
            # 'learning_rate': [i * 0.001 for i in range(1, 10)],
            # 'momentum': [i * 0.1 for i in range(1, 10)],
            # 'patience': range(1, 10),
            # 'min_improvement': [i * 0.001 for i in range(1, 10)],
            # 'validate_every': range(1, 10),
            # 'batch_size': range(1, 10),
            # 'hidden_dropout': [i * 0.1 for i in range(0, 10)],
            # 'input_dropout': [i * 0.1 for i in range(0, 10)]
            }

    classifier = ClassifierGeneticAlgorithm(pop_size=int(argv[1]),
                                            pdb_data=pdb_data,
                                            split=float(argv[3]),
                                            args=args,
                                            retain=float(argv[4]),
                                            random_select=float(argv[5]),
                                            mutate=float(argv[6]),
                                            verbose=True)

    classifier.run()


if __name__ == '__main__':
    main(sys.argv)
