'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import itertools
import sys
import climate

from synbiochem.utils import structure_utils as struct_utils
import holygrail
import synbiochem.ann


def regress(sample_size, struct_pattern, split, hidden_layers):
    '''Regression of phi/psi angles of peptides, specified by structure
    pattern as regexps.'''
    climate.enable_default_logging()

    # Get random peptides that match structure patterns from PDB:
    pdb_data = holygrail.get_pdb_data(sample_size, [struct_pattern], True)

    # Convert peptides to inputs, based on amino acid properties:
    x_data = holygrail.get_input_data([i[0] for v in pdb_data.values()
                                       for i in v])

    # Convert classess to outputs (these are based on structure patterns):
    y_data = [_get_phi_psi_data(v[2][0], v[2][1], v[3])
              for sublist in pdb_data.values()
              for v in sublist]

    # Randomise input and output data order:
    x_data, y_data = synbiochem.ann.randomise_order(x_data, y_data)

    # Split data into training and classifying:
    ind = int(split * len(x_data))

    # Perform regression:
    regressor = synbiochem.ann.Regressor()
    regressor.train(x_data[:ind], y_data[:ind], hidden_layers=[hidden_layers])

    y_pred = regressor.predict(x_data[ind:])

    return y_data[ind:], y_pred


def _get_proximity_data(pdb_id):
    '''Gets proximity_data for deep learning.'''
    all_sequences = struct_utils.get_sequences(pdb_id)
    all_proximity_data = struct_utils.calc_proximities(pdb_id)

    # Ensure data consistency in terms of data lengths:
    assert len(all_sequences) == len(all_proximity_data)
    assert all([len(all_sequences[idx]) == len(all_proximity_data[idx])
                for idx in range(len(all_sequences))])

    prox_output = [list(itertools.chain.from_iterable(proximities))
                   for proximities in all_proximity_data]

    return all_sequences, prox_output


def _get_phi_psi_data(pdb_id, chain, subrange):
    '''Gets input/output for deep learning.'''
    all_phi_psi_data = struct_utils.get_phi_psi_data(pdb_id, chain)

    phi_psi_output = [list(itertools.chain.from_iterable(lst))
                      for lst in all_phi_psi_data]

    return phi_psi_output[0][slice(*[x * 2 for x in subrange])]


def _get_sub_square_matrix(idx, lngth, size):
    '''Returns square submatrix of length lngth from square matrix of given
    size.'''
    return [((idx + r) * size) + idx + c for r in range(lngth)
            for c in range(lngth)]


def main(argv):
    '''main method.'''
    regression = regress(int(argv[1]),
                         argv[4],
                         float(argv[2]),
                         int(argv[3]))

    for output in regression:
        print output


if __name__ == '__main__':
    main(sys.argv)
