'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import itertools
import random

from synbiochem.utils import structure_utils as struct_utils
import holygrail


def get_regression_data(max_ids, num_samples, nmer_len):
    '''Gets random data for regression analyses.'''
    for pdb_id in struct_utils.get_pdb_ids(max_ids):
        try:
            for data in _sample_regression_data(pdb_id, num_samples,
                                                nmer_len):
                print '\t'.join([str(datum) for datum in data])
        except AssertionError, err:
            print err


def _sample_regression_data(pdb_id, num_samples, nmer_len):
    '''Samples learning data for deep learning.'''
    pdb_id, all_seqs, inpt, prox_output, phi_psi_output = \
        _get_regression_data(pdb_id)

    results = []

    num_aa_props = len(holygrail.AMINO_ACID_PROPS['A'])

    for _ in range(num_samples):
        chn = int(random.random() * len(all_seqs))
        start = int(random.random() * (len(all_seqs[chn]) - nmer_len))
        rnge = range(start, start + nmer_len)
        result = [pdb_id,
                  ''.join([all_seqs[chn][i] for i in rnge]),
                  inpt[chn][start * num_aa_props:(start + nmer_len) *
                            num_aa_props],
                  [prox_output[chn][i]
                   for i in _get_sub_square_matrix(start, nmer_len,
                                                   len(all_seqs[chn]))],
                  phi_psi_output[chn][start * 2:(start + nmer_len) * 2]]

        results.append(result)

    return results


def _get_regression_data(pdb_id):
    '''Gets input/output for deep learning.'''
    all_sequences = struct_utils.get_sequences(pdb_id)
    all_proximities = struct_utils.calc_proximities(pdb_id)
    all_phi_psi_data = struct_utils.get_phi_psi_data(pdb_id)

    # Ensure data consistency in terms of data lengths:
    assert len(all_sequences) == len(all_proximities) == \
        len(all_phi_psi_data)
    assert all([len(all_sequences[idx]) == len(all_proximities[idx]) ==
                len(all_phi_psi_data[idx])
                for idx in range(len(all_sequences))])

    prox_output = [list(itertools.chain.from_iterable(proximities))
                   for proximities in all_proximities]

    phi_psi_output = [list(itertools.chain.from_iterable(lst))
                      for lst in all_phi_psi_data]

    return pdb_id, all_sequences, holygrail.get_input_data(all_sequences), \
        prox_output, phi_psi_output


def _get_sub_square_matrix(idx, lngth, size):
    '''Returns square submatrix of length lngth from square matrix of given
    size.'''
    return [((idx + r) * size) + idx + c for r in range(lngth)
            for c in range(lngth)]
