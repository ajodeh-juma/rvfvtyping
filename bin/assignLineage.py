#!/usr/bin/env python

import argparse
from lineageFinder import LineageFinder, prepare_tree



def main(args=None):
    parser = argparse.ArgumentParser(
        description="Searches a tree for a taxon and returns that taxon's lineage and bootstrap label for the "
                    "appropriate node in csv format. All other taxa should be labeled with lineage in the taxon "
                    "name. The location of the lineage in the name is specified with the separator and index "
                    "arguments.",
        usage="--separator \"|\" --index 2 --taxon my|fav|taxon -i my.tree -o my.csv --collapse_to_polytomies",
    )

    required = parser.add_argument_group("Required")
    required.add_argument(
        "-i",
        "--input",
        metavar='input.tree',
        dest="input",
        type=str,
        required=True,
        help='The input tree file. Format can be specified with the format flag.')
    required.add_argument(
        "-o",
        "--output",
        metavar='output.*',
        dest="output",
        type=str,
        required=True,
        help='The output file')

    required.add_argument(
        "--index",
        dest="index",
        type=int,
        required=True,
        help="The index of the trait to reconstruct when the tip label is split by the  separator"
    )

    required.add_argument(
        "-s",
        "--separator",
        dest="separator",
        type=str,
        required=True,
        help="optional separator used to get trait"
    )

    required.add_argument(
        "-t",
        "--taxon",
        dest='taxon',
        type=str,
        required=True,
        help='The tip label to get')

    parser.add_argument(
        '--format',
        dest='format',
        action='store',
        default="nexus",
        choices=['nexus', 'newick', 'nexml'],
        help='what format is the tree file. This is passed to dendropy. default is \'nexus\'')

    parser.add_argument(
        "-c",
        "--collapse_to_polytomies",
        dest='collapse',
        type=float,
        help='A optional flag to collapse branches with length <= the input before running the rule.')

    args = parser.parse_args()
    tree = prepare_tree(args)
    finder = LineageFinder(tree, args.taxon, args.index, args.separator)
    with open(args.output, "w") as output_file:
        output_file.write("%s,%s,%s,%s" % (args.taxon, finder.lineage, str(finder.alrt), str(finder.boot)))


if __name__ == "__main__":
    main()