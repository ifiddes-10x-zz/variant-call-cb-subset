"""
Parses the inputs provided by cellranger-dna and walks the tree, producing barcode sets as input for FILTER_BAM

"""
import martian
import pandas as pd
from scipy.cluster import hierarchy
from crdna.utils import load_h5


__MRO__ = """
stage PARSE_INPUTS(
    in csv barcodes,
    in h5 tree_data,
    in csv barnyard,
    in int min_cells,
    out string[] barcode_subsets
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
    barcode_map = {int(x.rsplit('_', 1)[-1]): y for x, y in zip(cell_df.cell_id, cell_df.BC)}

    # load the set of barcodes to be investigated
    if args.barcodes is not None:
        relevant_barcodes = {x.rstrip() for x in open(args.barcodes)}
        subset_df = cell_df[cell_df.BC.isin(relevant_barcodes)]
        assert len(subset_df) > 0
        if len(relevant_barcodes) != len(subset_df):
            missing = relevant_barcodes - set(subset_df.BC)
            martian.log_warn('Cell barcodes in whitelist are not cell barcodes:\n{}'.format(', '.join(missing)))
        relevant_cell_ids = {int(x.rsplit('_', 1)[-1]) for x in subset_df.cell_id}
    else:
        # all cells
        relevant_barcodes = None
        relevant_cell_ids = set(range(len(cell_df)))

    # load the linkage matrix
    Z = load_h5(args.tree_data, "Z").values

    # convert the linkage matrix into a ClusterNode structure
    rootnode, nodelist = hierarchy.to_tree(Z, rd=True)

    # find the root node of our given barcode subset, if it is not None
    if relevant_barcodes is not None:
        rootnode = find_root_node(nodelist, relevant_cell_ids)

    # now find all barcode subsets
    barcode_subsets = []
    for n in find_all_internal_nodes(rootnode):
        bcs = {barcode_map[x] for x in n.pre_order()}
        if len(bcs) >= args.min_cells:
            barcode_subsets.append(','.join(bcs))
    assert len(barcode_subsets) > 0
    outs.barcode_subsets = barcode_subsets


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


def find_root_node(nodelist, relevant_cell_ids):
    """
    Find the lowest common ancestor of relevant_indices
    """
    children_sets = {node: set(node.pre_order()) for node in nodelist}
    relevant_children_sets = {node: s for node, s in children_sets.iteritems()
                              if len(s & relevant_cell_ids) == len(relevant_cell_ids)}
    return sorted(relevant_children_sets.iteritems(), key=lambda x: len)[0][0]