from collections import defaultdict
import re

import argh

from helpers import load_exons


def main(input_filename, organism_id, output_filename):
    exons_by_gene_id = load_exons(input_filename, organism=organism_id)
    # now, we possibly have exons from multiple transcripts under the same gene_id
    # there are two formats in which the exons are named: '<gene_id>-R<transcript_id>-E<exon_id>' or 'exon_<something>_t<transcript_id>-E<exon_id>'
    # we want to keep only the longest transcript for each gene_id. Longest = the sum of lengths of the exons is the longest
    exon_format_1 = re.compile(
        r"^(.+)-R(\w+)-E(\w+)$"
    )  # <gene_id>-R<transcript_id>-E<exon_id>'
    exon_format_2 = re.compile(
        r"^exon_(.+)_t(\w+)-E(\w+)$"
    )  # exon_<something>_t<transcript_id>-E<exon_id>

    filtered_exons = []

    exons_filtered_by_gene_id = {}
    for gene_id, exons in exons_by_gene_id.items():
        # let's split the exons by the transcript
        transcripts = defaultdict(list)
        for exon in exons:
            exon_full_id = exon.exon_id

            if (m := exon_format_1.match(exon_full_id)) is not None:
                transcript_id = m.group(2)
                exon_num = m.group(3)
            elif (m := exon_format_2.match(exon_full_id)) is not None:
                transcript_id = m.group(2)
                exon_num = m.group(3)
            else:
                raise ValueError(f"Unknown exon format: {exon_full_id}")
            transcripts[transcript_id].append(
                exon._replace(exon_id=f"{gene_id}-R{transcript_id}-E{exon_num}")
            )
        # let's compute the combined length of the exons for each transcript
        transcript_lengths = {}
        for transcript_id, exons in transcripts.items():
            transcript_lengths[transcript_id] = sum(
                exon.end - exon.start for exon in exons
            )
        # let's find the longest transcript
        longest_transcript_id = max(transcript_lengths, key=transcript_lengths.get)
        filtered_exons.extend(transcripts[longest_transcript_id])

    filtered_exons = sorted(filtered_exons)

    with open(output_filename, "w") as output_file:
        for exon in filtered_exons:
            output_file.write(
                f"{exon.chromosome}\t{exon.start}\t{exon.end}\t{exon.exon_id}\t{0}\t{exon.strand}\n"
            )


if __name__ == "__main__":
    argh.dispatch_command(main)
