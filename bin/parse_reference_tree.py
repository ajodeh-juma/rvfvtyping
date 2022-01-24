#!/usr/bin/env python

import dendropy
import warnings
from collections import Counter
from ete3 import Tree


def parse_tree(filename, schema, threshold=80):
    """

    :param filename:
    :param schema:
    :param threshold:
    :return:
    """

    tree = dendropy.Tree.get(path=filename, schema=schema.lower(), preserve_underscores=True)
    # Remove nodes if their support value or branch length are below the given threshold
    for node in tree.postorder_node_iter():
        if node.label is None:
            continue
        else:
            value = node.label.split("/")[1]
            print(node, value, node.label)
    #     try:
    #         value = node.label.split("/")[1]
    #         if float(value) < threshold:
    #             for child in node.child_nodes():
    #                 child.parent_node = node.parent_node
    #             if node.parent_node:  # Ensure not working from the root
    #                 node.parent_node.remove_child(node)
    #                 tree.encode_bipartitions()
    #     except TypeError:
    #         continue
    # multifurcating = lambda x: True if len(x.child_nodes()) > 1 else False
    for nd in tree.preorder_node_iter():
        print(nd.description())
    #
    # for taxon in tree.taxon_namespace:
    #     print(taxon)
    t = Tree(filename, format=1)
    t.show()
    return tree


def parse_beast_tree(filename, schema):
    """

    :param filename:
    :param schema:
    :param threshold:
    :return:
    """

    tree = dendropy.Tree.get_from_path(src=filename, schema=schema.lower(), preserve_underscores=False)
    # Remove nodes if their support value or branch length are below the given threshold

    for node in tree.postorder_node_iter():
        print(node, node.taxon, node.posterior)
    return tree


if __name__ == '__main__':
    tree_file = "/Users/jjuma/Documents/PhD_RVF2019/pipelines/rvfvtyping/rvfv-typing/test/validation-rep-outdir/iqtree/repseq.align.fasta.treefile"
    tree = parse_tree(filename=tree_file, schema='newick', threshold=80)

    # tree_file = "/Users/jjuma/Documents/PhD_RVF2019/pipelines/rvfvtyping/rvfv-lineages/G2-glycoprotein_Lineages.align.dedup.MCC.tree"
    # tree = parse_beast_tree(filename=tree_file, schema="nexus")
