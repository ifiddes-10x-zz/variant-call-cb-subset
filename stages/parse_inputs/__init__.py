"""
Parses the inputs provided by cellranger-dna and walks the tree, producing barcode sets as input for FILTER_BAM

"""
import martian
import pandas as pd
from scipy.cluster import hierarchy
from crdna.utils import load_h5
from crdna.clustering import *


__MRO__ = """
stage PARSE_INPUTS(
    in csv barcodes,
    in h5 tree_data,
    in csv barnyard,
    in int min_cells,
    out string[] barcode_subsets,
    out string[] node_ids,
    src py "stages/parse_inputs",
)
"""

def main(args, outs):
    """
    Performs the splitting operation.
    """
    # load the barnyard dataframe -- this maps cell IDs into barcodes
    barnyard_df = pd.read_csv(args.barnyard)
    # extract relevant entries
    cell_df = barnyard_df[barnyard_df['cell_id'] != 'None']
    barcodes = np.array(cell_df['BC'])

    # load the linkage matrix
    Z = load_h5(args.tree_data, "Z").values
    desc = get_descendant_matrix(Z)

    # convert the linkage matrix into a ClusterNode structure
    rootnode, nodelist = hierarchy.to_tree(Z, rd=True)

    # load the set of barcodes to be investigated
    if args.barcodes is not None:
        relevant_barcodes = {x.rstrip() for x in open(args.barcodes)}
        subset_df = cell_df[cell_df.BC.isin(relevant_barcodes)]
        assert len(subset_df) > 0
        if len(relevant_barcodes) != len(subset_df):
            missing = relevant_barcodes - set(subset_df.BC)
            martian.log_warn('Cell barcodes in whitelist are not cell barcodes:\n{}'.format(', '.join(missing)))
        relevant_cell_ids = {int(x.rsplit('_', 1)[-1]) for x in subset_df.cell_id}
        relevant_cell_ids = np.array([1 if x in relevant_cell_ids else 0 for x in xrange(desc.shape[0])])
    else:
        # all cells
        relevant_barcodes = None
        relevant_cell_ids = np.array([1] * desc.shape[0])

    # find the root node of our given barcode subset, if it is not None
    if relevant_barcodes is not None:
        rootnode = nodelist[find_root_node(desc, relevant_cell_ids)]

    # now find all barcode subsets
    barcode_subsets = [','.join(barcodes[desc[rootnode.id, :]])]
    # record the node ID for the final VCF
    node_ids = [str(rootnode.id)]
    if args.min_cells is not None:
        for n in find_all_internal_nodes(rootnode):
            bcs = barcodes[desc[n.id, :]]
            if len(bcs) >= args.min_cells:
                barcode_subsets.append(','.join(bcs))
                node_ids.append(str(n.id))

    assert len(barcode_subsets) > 0
    outs.barcode_subsets = barcode_subsets
    outs.node_ids = node_ids


def find_all_internal_nodes(root):
    """
    Binary search, yields only internal nodes
    """
    left = root.get_left()
    if left is not None:
        if not left.is_leaf():
            yield left
        for n in find_all_internal_nodes(left):
            yield n
    right = root.get_right()
    if right is not None:
        if not right.is_leaf():
            yield right
        for n in find_all_internal_nodes(right):
            yield n


def find_root_node(desc, labels):
    ncells = desc.shape[1]
    node_ncell = desc[:, :ncells].sum(axis=1)
    node_ncell = sorted(enumerate(node_ncell), key=lambda x: x[1], reverse=True)
    root = -1
    for node, n in node_ncell:
        n_label = np.sum(labels[desc[node, :].nonzero()[0]])
        if n_label == n:
            root = node
            break
    return root if root >= 0 else None
