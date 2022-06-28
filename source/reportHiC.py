import strawC, numpy as np, random, argparse, gzip, pandas as pd


def readStrawCResult(result, bin_offset=0):
    i, j, v = [], [], []
    for r in result:
        i.append(r.binX)
        j.append(r.binY)
        v.append(r.counts)
    return np.array(i, dtype=int), np.array(j, dtype=int), np.array(v, dtype=float)


def main():

    parser = argparse.ArgumentParser(prog="base_maker", description="Extract Hi-C data")
    parser.add_argument("-f", "--hic-filename", help="path (local or remote) to Hi-C data in .hic format", type=str)
    parser.add_argument("-r", "--res", help="resolution to compute correlation at", type=int, default=1000)
    parser.add_argument("-m", "--max", help="maximum pairise distance", type=int, default=100000)
    parser.add_argument("-c", "--chrom", help="chromosome to process", type=str)
    parser.add_argument("-n", "--norm", type=str)
    parser.add_argument("-l", "--chrom-size-filename", type=str)
    parser.add_argument("-o", "--output-filename", help="path to a gzipped output file", type=str)
    args = parser.parse_args()

    # Read chromosome length
    c = args.chrom.replace("chr", "")
    csf = pd.read_csv(args.chrom_size_filename, sep="\t", header=None, names=["chrom", "length"])
    l = int(csf[csf.chrom == args.chrom]["length"])
    m = int(np.floor(l / 2 / args.res) * args.res)
    v = int(args.max / 2)
    print(l, m, v)

    with gzip.open(args.output_filename, "wt") as fout:
        # Read interactions in first half and then second half of the chromosome
        # (expensive and uneccessary to read interactions across the entire chromosome all at once)
        for i in range(2):
            coord = "%s:%d:%d" % (c, 0 if i == 0 else m - v, m + v if i == 0 else l)
            result = strawC.strawC("observed", args.norm, args.hic_filename, coord, coord, "BP", args.res)
            ia, ja, va = readStrawCResult(result, args.res)

            # Output to file
            for j in range(len(ia)):
                d = ja[j] - ia[j]
                if d >= args.res and d <= args.max:
                    fout.write("%s-%d-%d\t%.6f\n" % (args.chrom, ia[j], ja[j], va[j]))


main()
