import os, sys, argparse, random

import numpy as np

def parse_args():
    """ Return dictionary of command line arguments """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        usage=argparse.SUPPRESS)
    parser.add_argument('--catalog', type=str, dest='catalog', required=True,
                        help="""Path to snp catalog file""")
    parser.add_argument('--name', type=str, dest='name', required=True,
                        help="""A parameter specifying the name of the catalog""")
    parser.add_argument('--coords-dir', type=str, dest='coords_dir', required=True,
                        help="""A parameter specifying a path to the directory containing all coords files""")
    parser.add_argument('--out', type=str, dest='out', default='/dev/stdout',
                        help="""Path to output file""")

    return vars(parser.parse_args())


def read_catalog(fpath):
    snps = []
    samples = []
    catalog_vec = []

    count = 1
    with open(fpath) as fh:
        samples = next(fh).rstrip().split('\t')[1:]
        for line in fh:
            items = line.rstrip().split('\t')
            snps.append(items[0])
            count += 1
            for i in items[1:]:
                catalog_vec.append(np.uint8(i))

    sample_dict = dict()
    for i, sample in enumerate(samples):
        sample_dict[sample] = i

    return snps, sample_dict, np.reshape(catalog_vec, (len(snps), len(samples)))


def get_coords_paths(coords_dir):
    coords_paths = []
    for f in os.listdir(coords_dir):
        if 'coords' in f:
            fpath = coords_dir.rstrip('/') + '/' + f
            coords_paths.append(fpath)
        else:
            sys.stderr.write("skip {}\n".format(f))

    return coords_paths


def read_coords_test(coords_paths):
    ref_allele_dict = dict()

    for coords_path in coords_paths:
        ref_tag = coords_path.split('/')[-1].split('-')[0]
        genome_tag = coords_path.split('/')[-1].split('-')[-1].split('.')[0]

        if ref_tag == genome_tag:
            continue

        with open(coords_path) as fh:
            for i in range(5):
                next(fh)

            for line in fh:
                items = line.rstrip().split()
                start = int(items[0])
                end = int(items[1])
                chrom = items[11]
                snp_loci = [chrom + "||" + str(i) for i in range(start, end + 1)]

            sys.stderr.write("done {}\n".format(coords_path))

    return ref_allele_dict


def read_coords(coords_paths):
    ref_allele_dict = dict()

    for coords_path in coords_paths:
        ref_tag = coords_path.split('/')[-1].split('-')[0]
        genome_tag = coords_path.split('/')[-1].split('-')[-1].split('.')[0]

        if ref_tag == genome_tag:
            continue

        with open(coords_path) as fh:
            for i in range(5):
                next(fh)

            for line in fh:
                items = line.rstrip().split()
                start = int(items[0])
                end = int(items[1])
                chrom = items[11]

                for i in range(start, end + 1):
                    snp_locus = chrom + "||" + str(i)
                    if snp_locus not in ref_allele_dict:
                        ref_allele_dict[snp_locus] = []
                    ref_allele_dict[snp_locus].append(genome_tag)

            sys.stderr.write("done {}\n".format(coords_path))

    return ref_allele_dict


def infuse_refs(ref_allele_dict, catalog, samples, snps, output_path):
    for i, snp in enumerate(snps):
        id_items = snp.split('||')
        snp_locus = id_items[0] + '||' + id_items[1]

        if snp_locus not in ref_allele_dict:
            pass
        else:
            sample_indices = [samples[sample] for sample in ref_allele_dict[snp_locus]]
            alt_mask = catalog[i, :] == 1
            catalog[i, sample_indices] = 2
            miss_mask = catalog[i, :] == 0
            catalog[i, miss_mask] = 255
            catalog[i, sample_indices] = 0
            catalog[i, alt_mask] = 1

    sample_list = [samp_pair[0] for samp_pair in sorted(samples.items(), key=lambda kv: (kv[1], kv[0]))]

    with open(output_path, 'w') as fh:
        first_row = '\t'.join(['snp_id'] + sample_list)
        fh.write("{}\n".format(first_row))
        for i, row in enumerate(catalog):
            fh.write("{}\t{}\n".format(snps[i], '\t'.join(map(str, row))))


def main():
    args = parse_args()
    catalog_path = args["catalog"]
    coords_dir = args["coords_dir"]
    name = args["name"]
    output_path = args["out"]

    snps, samples, catalog = read_catalog(catalog_path)
    print(samples)

    coords_paths = get_coords_paths(coords_dir)
    ref_allele_dict = read_coords(coords_paths)

    infuse_refs(ref_allele_dict, catalog, samples, snps, output_path)


if __name__ == "__main__":
    main()
