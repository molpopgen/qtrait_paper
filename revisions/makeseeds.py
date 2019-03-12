import argparse
import sys
import numpy as np


def make_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument("num_seed", type=int)
    parser.add_argument("outfile", type=str)
    parser.add_argument("seed", type=int)

    return parser


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    np.random.seed(args.seed)
    seeds = []
    for i in range(args.num_seed):
        candidate = np.random.randint(
            np.iinfo(np.uint32).max, dtype=np.uint32)
        while candidate in seeds:
            candidate = np.random.randint(
                np.iinfo(np.uint32).max, dtype=np.uint32)
        seeds.append(candidate)

    with open(args.outfile, 'w') as f:
        for i in seeds:
            f.write("{}\n".format(i))
