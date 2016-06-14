'''
HolyGrail (c) University of Manchester 2015

HolyGrail is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-arguments
import random

import synbiochem.utils.sequence_utils as seq_utils
import synbiochem.utils.structure_utils as struct_utils


def sample_seqs(sample_size, struct_sets, length, min_hamming=3):
    '''Sample sequence and structure data.'''
    seqs = {}
    hammings = {}

    for struct_set in struct_sets:
        seqs[struct_set] = []

        while len(seqs[struct_set]) < sample_size:
            matches_required = sample_size - len(seqs[struct_set])
            matches = random.sample(struct_utils.get_pep_structs(struct_set,
                                                                 length),
                                    matches_required)

            if min_hamming > 0:
                _filter_hammings(matches, sample_size, seqs, struct_set,
                                 min_hamming, hammings)
            else:
                seqs[struct_set] = matches

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
                hamming = seq_utils.get_hamming(seq, match)
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
