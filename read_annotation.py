import requests, sys

def get_rtf_data(rtf):
    exon_seqs = {}
    iso_per_gene = {}
    exons_per_iso = {}

    print('this is running')

    server = "https://rest.ensembl.org"

    table = open(rtf, 'r')
    lines = table.readlines()
    table.close()

    e = open('exons.fasta', 'w')
    a = open('annotation_data.csv', 'w')

    lc = 0

    for line in lines:
        print('Line ' + str(lc))
        lc += 1
        info = line[:-1].split("\t")
        type = info[2]
        tags = info[8][:-1].split("; ")
        tag_dict = {}
        for tag in tags:
            tag_info = tag.split(" ")
            tag_dict[tag_info[0]] = tag_info[1].replace('"', '')

        gene_id = tag_dict['gene_id'].split(".")[0]

        if type == 'exon':
            exon_id = tag_dict['exon_id'].split(".")[0]
            transcript_id = tag_dict['transcript_id'].split(".")[0]

            if gene_id not in iso_per_gene:
                iso_per_gene[gene_id] = set()
            if transcript_id not in exons_per_iso:
                exons_per_iso[transcript_id] = set()

            iso_per_gene[gene_id].add(transcript_id)
            exons_per_iso[transcript_id].add(exon_id)

            if exon_id not in exon_seqs:
                ext = "/sequence/id/" + exon_id + "?type=genomic"
                r = requests.get(server + ext, headers={"Content-Type": "text/plain"})

                if not r.ok:
                    r.raise_for_status()
                    sys.exit()

                exon_seqs[exon_id] = r.text

    for exon in exon_seqs:
        e.write(">" + exon + "\n")
        e.write(exon_seqs[exon] + "\n")

    a.write('gene,transcript,exon\n')
    for gene_id in iso_per_gene:
        for iso in iso_per_gene[gene_id]:
            for exon in exons_per_iso[iso]:
                a.write(gene_id + "," + iso + "," + exon + "\n")

    a.close()
    e.close()

    return

def main():
    get_rtf_data('liqa_example.gtf')

if __name__ == "__main__":
    main()
