'''
HolyGrail (c) University of Manchester 2015

HolyGrail is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys

import climate

import holygrail
import synbiochem.ann


def classify(sample_size, struct_patterns, split, hidden_layers):
    '''Classification of peptides, specified by structure patterns as
    regexps.'''
    climate.enable_default_logging()

    # Get random peptides that match structure patterns from PDB:
    pdb_data = holygrail.get_pdb_data(sample_size, struct_patterns)

    # Convert peptides to inputs, based on amino acid properties:
    x_data = holygrail.get_input_data([i[0] for v in pdb_data.values()
                                       for i in v])

    # Convert classess to outputs (these are based on structure patterns):
    y_data = [i for k, v in pdb_data.iteritems() for i in [k] * len(v)]

    # Randomise input and output data order:
    x_data, y_data = synbiochem.ann.randomise_order(x_data, y_data)

    # Split data into training and classifying:
    ind = int(split * len(x_data))

    # Perform classification:
    classifier = synbiochem.ann.Classifier()
    classifier.train(x_data[:ind], y_data[:ind], hidden_layers=[hidden_layers])

    return classifier.classify(x_data[ind:], y_data[ind:])[1:]


def main(argv):
    '''main method.'''
    classification = classify(int(argv[1]),
                              argv[4:],
                              float(argv[2]),
                              int(argv[3]))

    for output in classification:
        print output


if __name__ == '__main__':
    main(sys.argv)
