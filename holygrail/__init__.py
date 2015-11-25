'''
HolyGrail (c) University of Manchester 2015

HolyGrail is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import itertools
import random
import sys

import climate

import synbiochem.ann


# KD Hydrophobicity, EIIP, Helix, Sheet, Turn
AA_PROPS = {
    'A': [0.66, 0.336, 0.870, 0.399, 0.326],
    'R': [0.1, 0.707, 0.597, 0.518, 0.569],
    'N': [0.189, 0.123, 0.391, 0.276, 0.9],
    'D': [0.189, 0.9, 0.645, 0.276, 0.855],
    'C': [0.722, 0.625, 0.554, 0.608, 0.719],
    'Q': [0.189, 0.582, 0.691, 0.561, 0.59],
    'E': [0.189, 0.137, 0.9, 0.1, 0.403],
    'G': [0.464, 0.132, 0.213, 0.466, 0.9],
    'H': [0.216, 0.253, 0.691, 0.434, 0.569],
    'I': [0.9, 0.1, 0.589, 0.870, 0.1],
    'L': [0.838, 0.1, 0.814, 0.668, 0.251],
    'K': [0.153, 0.335, 0.755, 0.345, 0.61],
    'M': [0.669, 0.621, 0.828, 0.568, 0.263],
    'F': [0.749, 0.699, 0.683, 0.703, 0.263],
    'P': [0.358, 0.225, 0.1, 0.234, 0.883],
    'S': [0.429, 0.625, 0.323, 0.518, 0.842],
    'T': [0.438, 0.696, 0.391, 0.756, 0.576],
    'W': [0.42, 0.447, 0.622, 0.708, 0.576],
    'Y': [0.384, 0.427, 0.335, 0.746, 0.691],
    'V': [0.873, 0.136, 0.572, 0.9, 0.141]
}


def get_input_data(all_sequences):
    '''Returns input data for machine-learning problems.'''
    return [list(itertools.chain.from_iterable([AA_PROPS[am_acid]
                                                if am_acid in AA_PROPS
                                                else [0.5] * len(AA_PROPS['A'])
                                                for am_acid in sequences]))
            for sequences in all_sequences]
