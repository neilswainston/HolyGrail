'''
HolyGrail (c) University of Manchester 2015

HolyGrail is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-few-public-methods
import sys

import holygrail.data
import synbiochem.ann


# import climate
class Classifier(object):
    '''Class to represent a Classifier of secondary structure.'''

    def __init__(self, pdb_data, split, scale=(0.1, 0.9)):
        '''Constructor.'''
        # climate.enable_default_logging()

        # Convert peptides to inputs, based on amino acid properties:
        x_data = holygrail.get_input_data([i[0] for v in pdb_data.values()
                                           for i in v],
                                          scale=scale)

        # Convert classes to outputs (these are based on structure patterns):
        y_data = [i for k, v in pdb_data.iteritems() for i in [k] * len(v)]

        # Randomise input and output data order:
        x_data, y_data = synbiochem.ann.randomise_order(x_data, y_data)

        # Split data into training and testing:
        ind = int(split * len(x_data))
        self.__x_data_train = x_data[:ind]
        self.__y_data_train = y_data[:ind]
        self.__x_data_test = x_data[ind:]
        self.__y_data_test = y_data[ind:]

    def classify(self, hidden_layers=None, input_noise=0.0, hidden_noise=0.0,
                 learning_rate=0.01, momentum=0.5, patience=5,
                 min_improvement=0.005, validate_every=1, batch_size=5,
                 hidden_dropout=0.0, input_dropout=0.0,
                 aa_props_filter=(2**holygrail.NUM_AA_PROPS - 1)):
        '''Classification of peptides, specified by structure patterns as
        regexps.'''

        x_data_train = _filter_x_data(self.__x_data_train, aa_props_filter)
        x_data_test = _filter_x_data(self.__x_data_test, aa_props_filter)

        # Perform classification:
        classifier = synbiochem.ann.Classifier(x_data_train,
                                               self.__y_data_train)
        classifier.train(hidden_layers=hidden_layers, input_noise=input_noise,
                         hidden_noise=hidden_noise,
                         learning_rate=learning_rate, momentum=momentum,
                         patience=patience, min_improvement=min_improvement,
                         validate_every=validate_every, batch_size=batch_size,
                         hidden_dropout=hidden_dropout,
                         input_dropout=input_dropout)

        return classifier.classify(x_data_test, self.__y_data_test)

    def get_y_data_test(self):
        '''Returns the test output data.'''
        return self.__y_data_test


def _filter_x_data(x_data, aa_props_filter):
    '''Filter x_data (effectively select amino acid parameters).'''
    expand_binary = list(reversed([int(d)
                                   for d in str(bin(aa_props_filter))[2:]]))
    flt = [i for i, x in enumerate(expand_binary) if x == 1]
    return [[v for i, v in enumerate(values)
             if i % holygrail.NUM_AA_PROPS in flt]
            for values in x_data]


def main(argv):
    '''main method.'''
    # Get random peptides that match structure patterns from PDB:
    pdb_data, hammings = holygrail.data.sample_seqs(int(argv[1]),
                                                    argv[4:])

    print hammings

    classifier = Classifier(pdb_data, float(argv[2]))
    classification = classifier.classify(hidden_layers=[int(argv[3])])

    for test, pred in zip(classifier.get_y_data_test(), classification[0]):
        print '\t'.join([test, pred, str(test == pred)])

    for output in classification[1:]:
        print output


if __name__ == '__main__':
    main(sys.argv)
