#! /usr/bin/env python
# Python 2/3 compatibility
import sys
import math
import pip
import os
try:
    import more_itertools
except ImportError:
    os.system('pip install more_itertools')


def binary_search(start, end, informative_sites):
    matches = []
    query_start = 0
    query_end = len(informative_sites) - 1
    query_start_prev = -1
    query_end_prev = -1
    while len(matches) <= 0 and query_end > -1:
        # if the query region converges, no sites in read
        if query_start > query_end:
            break
        # if the query region isn't changing, no sites in read
        if (query_start == query_start_prev) and (query_end == query_end_prev):
            break

        query_start_prev = query_start
        query_end_prev = query_end
        query_pos = int((query_end + query_start) / 2)
        # if the query position is between the start and the end of the read, we've arrived
        if start <= informative_sites[query_pos]["pos"] < end:
            # keep the one that we found
            matches.append(informative_sites[query_pos])

            # check the next ones to see if they also fit the bill
            for isite in informative_sites[query_pos + 1:]:
                if start <= isite["pos"] <= end:
                    matches.append(isite)
                else:
                    break
            # and the previous ones
            for isite in informative_sites[:query_pos][::-1]:
                if start <= isite["pos"] <= end:
                    matches.append(isite)
                else:
                    break
            break
        elif informative_sites[query_pos]["pos"] > start:
            # move left, position is too high
            query_end = query_pos - 1
        elif informative_sites[query_pos]["pos"] < start:
            # move right, position is too low
            query_start = query_pos + 1
    return matches


def match_informative_sites(reads, informative_sites):
    """
    Given a list of pysam reads,
    and a dictionary of lists of informative sites,
    return the informative sites that are in the reads.
    Pasing in chrom because the pysam object keeps getting it wrong
    """
    matches = {}

    # reads are stored by whether they support ref or alt alleles
    for ref_alt in reads:
        matches[ref_alt] = []

        for read in reads[ref_alt]:
            site_matches = binary_search(
                read.reference_start, read.reference_end, informative_sites
            )

            if len(site_matches) > 0:
                ref_parents = set()
                alt_parents = set()
                for match in site_matches:
                    ref_parents.add(match["ref_parent"])
                    alt_parents.add(match["alt_parent"])
                if len(ref_parents) != 1 or len(alt_parents) != 1:
                    continue  # indicates a read that supports both parents, which doesn't work biologically
                match_info = {"matches": site_matches, "read": read}
                matches[ref_alt].append(match_info)
    return matches


def verify(region_chrom, dad, mom, bamfile, reads, informative_sites):
    """
    The idea is to loop over all informative sites first, and then loop over the reads for each informative site.
    If reads at this informative site can prove which parent the denovo mutation came from, then record it. The idea is
    basically the same as the naked eye to find the denovo mutation
    :param region_chrom: Chromosome name
    :param dad: dad
    :param mom: mom
    :param bamfile: bam file
    :param reads: All reads detected by the current denovo mutation site sequencing
    :param informative_sites: All informative sites for the current denovo mutation
    :return: Returns the content of each line of output
    """
    site_support = {}
    for site in informative_sites:
        try:
            bam_iter = bamfile.fetch(region_chrom, site["pos"], site["pos"] + 1)
        except ValueError:
            chrom = (
                region_chrom.strip("chr")
                if "chr" in region_chrom
                else "chr" + region_chrom
            )
            bam_iter = bamfile.fetch(chrom, site["pos"], site["pos"] + 1)
        site_depth = more_itertools.ilen(bam_iter)
        seq_error = math.ceil(site_depth * 0.15)
        denovo_muta = {}
        for ref_alt in reads:
            for read in reads[ref_alt]:
                try:
                    infosite_pos = read.get_reference_positions(full_length=True).index(int(site["pos"]))
                except ValueError:
                    continue
                kid_allele = read.query_sequence[infosite_pos]

                read_origin = ""
                if kid_allele == site["ref_allele"]:
                    read_origin = "ref_parent"
                elif kid_allele == site["alt_allele"]:
                    read_origin = "alt_parent"
                else:
                    continue

                if read_origin == "ref_parent":
                    if ref_alt == "ref":
                        if site["alt_parent"] not in denovo_muta:
                            denovo_muta[site["alt_parent"]] = []
                        denovo_muta[site["alt_parent"]].append(read.query_name)
                    else:
                        if site["ref_parent"] not in denovo_muta:
                            denovo_muta[site["ref_parent"]] = []
                        denovo_muta[site["ref_parent"]].append(read.query_name)
                else:
                    if ref_alt == "ref":
                        if site["ref_parent"] not in denovo_muta:
                            denovo_muta[site["ref_parent"]] = []
                        denovo_muta[site["ref_parent"]].append(read.query_name)
                    else:
                        if site["alt_parent"] not in denovo_muta:
                            denovo_muta[site["alt_parent"]] = []
                        denovo_muta[site["alt_parent"]].append(read.query_name)

        if len(denovo_muta.keys()) == 2:
            if len(denovo_muta[dad]) > seq_error and len(denovo_muta[dad]) > len(denovo_muta[mom]):
                # The third-generation sequencing error rate is 15%. If there are enough reads to prove that
                # the denovo mutation originated from either parent, and the number of reads is greater than the number
                # of sequencing errors at the current locus, it can be proved that this locus can support
                # the denovo mutation.
                denovo_muta["denovo_parent"] = dad
            elif len(denovo_muta[mom]) > seq_error and len(denovo_muta[mom]) > len(denovo_muta[dad]):
                denovo_muta["denovo_parent"] = mom
            else:
                continue
        elif len(denovo_muta.keys()) == 1:
            denovo_muta["denovo_parent"] = list(denovo_muta.keys()).pop()
        else:
            continue
        site_support[site["pos"]] = denovo_muta

    denovo_result = [site_support[item]["denovo_parent"] for item in site_support.keys()]
    if len(denovo_result) == 0:
        return []
    dad_ratio = denovo_result.count(dad) / len(denovo_result)
    mom_ratio = denovo_result.count(mom) / len(denovo_result)

    if dad_ratio > 0.8:
        # 0.8 is 4:1, which means that we believe that four of the five information sites prove that
        # the denovo mutation originated from the father, which has a high enough accuracy and can be modified
        # artificially.
        half_output = [dad, mom, str(denovo_result.count(dad)), 'READBACKED']
        output_site = []
        output_read = []
        for key in site_support.keys():
            if site_support[key]["denovo_parent"] == dad:
                output_site.append(str(key))
                output_read.append(site_support[key][dad].pop())
        half_output.append(','.join(output_site))
        half_output.append(','.join(output_read))
        half_output.extend(['-', '-'])
        return half_output
    elif mom_ratio > 0.8:
        half_output = [mom, dad, str(denovo_result.count(mom)), 'READBACKED']
        output_site = []
        output_read = []
        for key in site_support.keys():
            if site_support[key]["denovo_parent"] == mom:
                output_site.append(str(key))
                output_read.append(site_support[key][mom].pop())
        half_output.append(','.join(output_site))
        half_output.append(','.join(output_read))
        half_output.extend(['-', '-'])
        return half_output
    else:
        half_output = ["ambiguous", "ambiguous", str(0), 'READBACKED']
        output_site_dad = []
        output_read_dad = []
        output_site_mom = []
        output_read_mom = []
        for key in site_support.keys():
            if site_support[key]["denovo_parent"] == dad:
                output_site_dad.append(str(key))
                output_read_dad.append(site_support[key][dad].pop())
            else:
                output_site_mom.append(str(key))
                output_read_mom.append(site_support[key][mom].pop())
        half_output.append(','.join(output_site_dad))
        half_output.append(','.join(output_read_dad))
        half_output.append(','.join(output_site_mom))
        half_output.append(','.join(output_read_mom))
        return half_output


if __name__ == "__main__":
    sys.exit("Import this as a module")
