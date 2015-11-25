'''
HolyGrail (c) University of Manchester 2015

HolyGrail is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys

import climate

from synbiochem.utils import structure_utils as struct_utils
import holygrail.ann


def get_classif_data(sample_size, struct_patterns):
    '''Gets random data for classification analyses.'''
    return {struct_pattern: struct_utils.sample_seqs(sample_size,
                                                     struct_pattern)
            for struct_pattern in struct_patterns}


def main(argv):
    '''main method.'''
    climate.enable_default_logging()

    classif_data = get_classif_data(int(argv[1]), argv[3:])

    x_data = holygrail.get_input_data([i for v in classif_data.values()
                                       for i in v])
    y_data = [i for k, v in classif_data.iteritems() for i in [k] * len(v)]

    x_data, y_data = holygrail.ann.randomise_order(x_data, y_data)

    # Split data into training and classifying:
    ind = int(0.8 * len(x_data))

    classifier = holygrail.ann.Classifier()
    classifier.train(x_data[:ind], y_data[:ind], hidden_layers=[int(argv[2])])

    for output in classifier.classify(x_data[ind:], y_data[ind:])[1:]:
        print output

if __name__ == '__main__':
    main(sys.argv)
