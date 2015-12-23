'''
HolyGrail (c) University of Manchester 2015

HolyGrail is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import matplotlib.pyplot
import sys

import synbiochem.utils.structure_utils as struct_utils


def main(argv):
    '''main method.'''
    seqs, hammings = struct_utils.sample_seqs(int(argv[1]),
                                              argv[4:],
                                              min_hamming=int(argv[2]))

    matplotlib.pyplot.bar(hammings.keys(), hammings.values())
    print 'Hamming distances: ' + str(hammings)
    print 'Number of Hamming distances: ' + \
        str(sum([v for v in hammings.values()]))

    with open(argv[3], 'w') as outfile:
        for struct_patt, values in seqs.iteritems():
            for seq in values:
                outfile.write('\t'.join([struct_patt] +
                                        [str(v) for v in seq]) + '\n')

if __name__ == '__main__':
    main(sys.argv)
