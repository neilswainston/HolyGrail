'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import itertools
import sys

from sklearn.metrics import mean_squared_error
import climate

from synbiochem.utils import sequence_utils, structure_utils
import holygrail
import synbiochem.ann


def regress(sample_size, struct_pattern, split, hidden_layers):
    '''Regression of phi/psi angles of peptides, specified by structure
    pattern as regexps.'''
    climate.enable_default_logging()

    x_data = []
    y_data = []

    while len(x_data) < sample_size:
        # Get random peptides that match structure patterns from PDB:
        pdb_data, _ = holygrail.data.sample_seqs(sample_size, [struct_pattern],
                                                 local_only=True)

        # Convert peptides to inputs, based on amino acid properties:
        curr_x_data = sequence_utils.get_aa_props([i[0]
                                                   for v in pdb_data.values()
                                                   for i in v])

        # Get torsion angles as outputs:
        curr_y_data = [_get_phi_psi_data(v[2][0], v[2][1], v[3])
                       for sublist in pdb_data.values()
                       for v in sublist]

        x_data.extend([d for i, d in enumerate(curr_x_data)
                       if len(curr_y_data[i]) / 2 == len(curr_y_data[i]) /
                       len(sequence_utils.AA_PROPS['A'])])

        y_data.extend([d for i, d in enumerate(curr_y_data)
                       if len(curr_y_data[i]) / 2 == len(curr_y_data[i]) /
                       len(sequence_utils.AA_PROPS['A'])])

    # Randomise input and output data order:
    x_data, y_data = synbiochem.ann.randomise_order(x_data[:sample_size],
                                                    y_data[:sample_size])

    return _run_regressor(split, x_data, y_data, hidden_layers)


def _get_proximity_data(pdb_id):
    '''Gets proximity_data for deep learning.'''
    all_sequences = structure_utils.get_sequences(pdb_id)
    all_proximity_data = structure_utils.calc_proximities(pdb_id)

    # Ensure data consistency in terms of data lengths:
    assert len(all_sequences) == len(all_proximity_data)
    assert all([len(all_sequences[idx]) == len(all_proximity_data[idx])
                for idx in range(len(all_sequences))])

    prox_output = [list(itertools.chain.from_iterable(proximities))
                   for proximities in all_proximity_data]

    return all_sequences, prox_output


def _get_phi_psi_data(pdb_id, chain, subrange):
    '''Gets input/output for deep learning.'''
    all_phi_psi_data = structure_utils.get_phi_psi_data(pdb_id, chain)

    phi_psi_output = [list(itertools.chain.from_iterable(lst))
                      for lst in all_phi_psi_data]

    return phi_psi_output[0][slice(*[x * 2 for x in subrange])]


def _get_sub_square_matrix(idx, lngth, size):
    '''Returns square submatrix of length lngth from square matrix of given
    size.'''
    return [((idx + r) * size) + idx + c for r in range(lngth)
            for c in range(lngth)]


def _run_regressor(split, x_data, y_data, hidden_layers):
    '''Runs the regressor job.'''
    # Split data into training and classifying:
    ind = int(split * len(x_data))

    # Perform regression:
    regressor = synbiochem.ann.Regressor(x_data[:ind], y_data[:ind])
    regressor.train(hidden_layers=[hidden_layers])
    y_pred = regressor.predict(x_data[ind:])

    print len(y_data[ind:])
    print len(y_pred)
    print [len(d) for d in y_data[ind:]]
    print [len(d) for d in y_pred]

    return mean_squared_error(y_data[ind:], y_pred)


def main(argv):
    '''main method.'''
    print regress(int(argv[1]),
                  argv[4],
                  float(argv[2]),
                  int(argv[3]))


if __name__ == '__main__':
    main(sys.argv)
