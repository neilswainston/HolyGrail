'''
HolyGrail (c) University of Manchester 2015

HolyGrail is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import random
import sys

import regex

import synbiochem.utils.sequence_utils as seq_utils
import synbiochem.utils.structure_utils as struct_utils


def sample_seqs(sample_size, struct_patts, min_hamming=3, local_only=False):
    '''Sample sequence and structure data.'''
    seqs = {}
    hammings = {}

    for struct_patt in struct_patts:
        seqs[struct_patt] = []

        while len(seqs[struct_patt]) < sample_size:
            matches = []
            pdb_ids = struct_utils.get_pdb_ids(
                sample_size, local_only=local_only)
            seq_struct = struct_utils.get_seq_struct(pdb_ids)

            for pdb_ids, seq_struct in seq_struct.iteritems():
                try:
                    for match in regex.finditer(struct_patt, seq_struct[1],
                                                overlapped=True):
                        matches.append([seq_struct[0][slice(*(match.span()))],
                                        seq_struct[1][slice(*(match.span()))],
                                        pdb_ids,
                                        match.span()])
                except TypeError:
                    '''Take no action, but accept that seq_struct[1]
                    (the secondary structure) is occasionally and inexplicably
                    None.'''

            _filter_hammings(matches, sample_size, seqs, struct_patt,
                             min_hamming, hammings)

    return seqs, hammings


def _filter_hammings(matches, sample_size, seqs, struct_patt, min_hamming,
                     hammings):
    '''Filters matches by Hamming distance.'''
    for match in random.sample(matches, min(len(matches), sample_size -
                                            len(seqs[struct_patt]))):
        add = True
        curr_hamms = []
        for vals in seqs.values():
            for seq in vals:
                hamming = seq_utils.get_hamming(seq[0], match[0])
                if hamming < min_hamming:
                    add = False
                    break
                else:
                    curr_hamms.append(hamming)

        if add:
            seqs[struct_patt].append(match)

            for hamming in curr_hamms:
                if hamming not in hammings:
                    hammings[hamming] = 0

                hammings[hamming] += 1


def main(argv):
    '''main method.'''
    seqs, hammings = sample_seqs(int(argv[1]),
                                 argv[4:],
                                 min_hamming=int(argv[2]))

    print 'Hamming distances: ' + str(hammings)
    print 'Total: ' + \
        str(sum([v for v in hammings.values()]))

    with open(argv[3], 'w') as outfile:
        for struct_patt, values in seqs.iteritems():
            for seq in values:
                outfile.write('\t'.join([struct_patt] +
                                        [str(v) for v in seq]) + '\n')

if __name__ == '__main__':
    main(sys.argv)
