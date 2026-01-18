from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import Pool, TimeoutError
from functools import partial
import time
import logging
import subprocess
import os
import shutil
import re
import json
from BCBio import GFF
from collections import OrderedDict
from Bio.Seq import translate


def cds_judgment(cds_sequence, parse_phase=True, keep_stop=False, return_cds=False):
    """make sure a cds seq is good for translate"""

    if parse_phase:
        phase = 0
        orf_dict = {}

        for i in range(3):
            cds_now = cds_sequence[i:]
            aa_seq = translate(cds_now, to_stop=False)
            if '*' in aa_seq:
                star_index = aa_seq.index('*')
                one_star_aa_seq = aa_seq[:star_index+1]
            else:
                one_star_aa_seq = aa_seq
            phase = i
            orf_dict[phase] = (phase, one_star_aa_seq,
                               cds_now[:len(one_star_aa_seq)*3])

        best_phase = sorted(orf_dict, key=lambda x: len(
            orf_dict[x][1]), reverse=True)[0]

        phase, one_star_aa_seq, cds_now = orf_dict[best_phase]

        if len(one_star_aa_seq) * 3 / len(cds_sequence) > 0.95:
            good_orf = True
        else:
            good_orf = False

        if not keep_stop and '*' in one_star_aa_seq:
            out_aa_seq = one_star_aa_seq[:-1]
            out_cds_now = cds_now[:-3]
        else:
            out_aa_seq = one_star_aa_seq
            out_cds_now = cds_now

        if good_orf and len(out_aa_seq) * 3 != len(out_cds_now):
            raise ValueError("cds length error")

        if return_cds:
            return good_orf, phase, out_aa_seq, out_cds_now
        else:
            return good_orf, phase, out_aa_seq
    else:
        aa_seq = translate(cds_sequence, to_stop=False)

        if '*' in aa_seq:
            star_index = aa_seq.index('*')
            one_star_aa_seq = aa_seq[:star_index+1]
        else:
            one_star_aa_seq = aa_seq

        good_orf = True if len(cds_sequence) % 3 == 0 and len(
            one_star_aa_seq) == len(cds_sequence) / 3 else False

        if not keep_stop and '*' in one_star_aa_seq:
            out_aa_seq = one_star_aa_seq[:-1]
            out_cds_now = cds_sequence[:-3]
        else:
            out_aa_seq = one_star_aa_seq
            out_cds_now = cds_sequence

        if good_orf and len(out_aa_seq) * 3 != len(out_cds_now):
            raise ValueError("cds length error")

        if return_cds:
            return good_orf, None, out_aa_seq, out_cds_now
        else:
            return good_orf, None, out_aa_seq


class ChrLoci(object):
    def __init__(self, chr_id=None, strand=None, start=None, end=None, sp_id=None):
        self.chr_id = chr_id
        self.sp_id = sp_id

        if strand is None:
            self.strand = strand
        elif strand == "+" or str(strand) == '1':
            self.strand = "+"
        elif strand == "-" or str(strand) == '-1':
            self.strand = "-"
        else:
            self.strand = None

        if end is not None and start is not None:
            self.start = min(int(start), int(end))
            self.end = max(int(start), int(end))
            self._range = (self.start, self.end)
            self.length = abs(self.end - self.start) + 1


class GenomeFeature(ChrLoci):
    def __init__(self, id=None, type=None, chr_loci=None, qualifiers={}, sub_features=None, chr_id=None, strand=None, start=None, end=None, sp_id=None):
        if chr_loci:
            super(GenomeFeature, self).__init__(chr_id=chr_loci.chr_id, strand=chr_loci.strand, start=chr_loci.start,
                                                end=chr_loci.end, sp_id=chr_loci.sp_id)
        else:
            super(GenomeFeature, self).__init__(
                chr_id=chr_id, strand=strand, start=start, end=end, sp_id=sp_id)

        self.id = id
        self.type = type
        if chr_loci:
            self.chr_loci = chr_loci
        else:
            self.chr_loci = ChrLoci(chr_id=chr_id, strand=strand,
                                    start=start, end=end, sp_id=sp_id)
        self.sub_features = sub_features
        self.qualifiers = qualifiers

    def sgf_len(self):
        self.sgf_len_dir = {i: 0 for i in list(
            set([sgf.type for sgf in self.sub_features]))}

        for sgf in self.sub_features:
            sgf_len = abs(sgf.start - sgf.end) + 1
            self.sgf_len_dir[sgf.type] += sgf_len

    def get_bottom_subfeatures(self):
        bgf_list = []
        if sgf is None or len(self.sub_features) == 0:
            bgf_list.append(sgf)
            return bgf_list
        else:
            for sgf in sgf.sub_features:
                bgf_list.extend(self.get_bottom_subfeatures(sgf))
        return bgf_list

    def __eq__(self, other):
        return self.id == other.id and self.chr_loci == other.chr_loci and self.type == other.type

    def __hash__(self):
        return hash(id(self))


def ft2cl(feature_location, chr_id):
    """
    create ChrLoci by FeatureLocation from BCBio
    """
    return ChrLoci(chr_id=chr_id, strand=feature_location.strand, start=feature_location.start + 1,
                   end=feature_location.end)


def sf2gf(sf, chr_id):
    """
    create GenomeFeature by SeqFeature from BCBio
    """
    sf_cl = ft2cl(sf.location, chr_id)
    gf = GenomeFeature(id=sf.id, type=sf.type, chr_loci=sf_cl)
    gf.qualifiers = sf.qualifiers

    # parse sub_feature
    if hasattr(sf, 'sub_features') and len(sf.sub_features) != 0:
        gf.sub_features = []
        for sub_sf in sf.sub_features:
            gf.sub_features.append(sf2gf(sub_sf, chr_id))

    return gf


def read_gff_file(gff_file):
    feature_dict = OrderedDict()

    no_id = 0
    with open(gff_file, 'r') as in_handle:
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                new_feature = sf2gf(feature, rec.id)
                if new_feature.type not in feature_dict:
                    feature_dict[new_feature.type] = OrderedDict()
                if new_feature.id == '':
                    new_feature.id = 'NoID_%d' % no_id
                    no_id += 1
                feature_dict[new_feature.type][new_feature.id] = new_feature

    return feature_dict


def section(inter_a, inter_b, int_flag=False, just_judgement=False):
    """
    get the section
    :param inter_a:
    :param inter_b:
    :return:
    """
    all = sorted(list(inter_a) + list(inter_b))
    deta = (all[1], all[2])
    if int_flag is False:
        if max(inter_a) >= min(inter_b) and max(inter_b) >= min(inter_a):
            If_inter = True  # Yes
        else:
            If_inter = False  # No

        if just_judgement:
            return If_inter
        else:
            return If_inter, deta

    else:
        if max(inter_a) - min(inter_b) >= -1 and max(inter_b) - min(inter_a) >= -1:
            If_inter = True  # Yes
        else:
            If_inter = False  # No

        if just_judgement:
            return If_inter
        else:
            return If_inter, deta


def merge_intervals(input_list, int=False):
    """
    a function that will merge overlapping intervals
    :param intervals: a list of tuples
                      e.g. intervals = [(1,5),(33,35),(40,33),(10,15),(13,18),(28,23),(70,80),(22,25),(38,50),(40,60)]
    :param int: if the data is all int
    :return: merged list
    """
    intervals = []
    for i in input_list:
        intervals.append(tuple(i))

    sorted_by_lower_bound = sorted(intervals, key=lambda tup: min(tup))
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if int is False:
                if min(higher) <= max(lower):
                    upper_bound = max(lower + higher)
                    # replace by merged interval
                    merged[-1] = (min(lower), upper_bound)
                else:
                    merged.append(higher)
            elif int is True:
                if min(higher) <= max(lower):
                    upper_bound = max(lower + higher)
                    # replace by merged interval
                    merged[-1] = (min(lower), upper_bound)
                elif max(lower) + 1 == min(higher):
                    upper_bound = max(lower + higher)
                    # replace by merged interval
                    merged[-1] = (min(lower), upper_bound)
                else:
                    merged.append(higher)
    return merged


def overturn(inter_list):
    """
    input a list of intervals and the function will overturn it to a list with gap of the intervals
    :param inter_list:
    :return: gap_list
    """
    inter_list = sorted(merge_intervals(inter_list, True))
    output_list = []
    last_right = 0
    for index in range(0, len(inter_list) + 1):
        if index == 0:
            output_list.append((float('-inf'), inter_list[index][0] - 1))
            last_right = inter_list[index][1]
        elif index == len(inter_list):
            output_list.append((last_right + 1, float('inf')))
        else:
            output_list.append((last_right + 1, inter_list[index][0] - 1))
            last_right = inter_list[index][1]
    return output_list


def interval_minus_set(target, bullets):
    if len(bullets) == 0:
        return [target]
    gaps = overturn(bullets)
    output_list = []
    for i in gaps:
        If_inter, deta = section(target, i)
        if If_inter:
            output_list.append(deta)
    return output_list


def mkdir(dir_name, keep=True):
    if keep is False:
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
        os.makedirs(dir_name)
    else:
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

    return dir_name


def rmdir(dir_name):
    if os.path.exists(dir_name):
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)
        else:
            os.remove(dir_name)


def cmd_run(cmd_string, cwd=None, retry_max=5, silence=True):
    if not silence:
        print("Running " + str(retry_max) + " " + cmd_string)
    p = subprocess.Popen(cmd_string, shell=True,
                         stderr=subprocess.PIPE, stdout=subprocess.PIPE, cwd=cwd)
    output, error = p.communicate()
    if not silence:
        print(error.decode())
    returncode = p.poll()
    if returncode == 1:
        if retry_max > 1:
            retry_max = retry_max - 1
            cmd_run(cmd_string, cwd=cwd, retry_max=retry_max)

    output = output.decode()
    error = error.decode()

    return (not returncode, output, error)


def read_fasta(file_name):
    seqdict = {}

    f = open(file_name, 'r')
    all_text = f.read()
    # info = string.split(all_text, '>') python2
    info = all_text.split('\n>')
    while '' in info:
        info.remove('')
    for i in info:
        # seq = string.split(i, '\n', 1) python2
        seq = i.split('\n', 1)
        seq[1] = re.sub(r'\n', '', seq[1])
        seq[1] = re.sub(r' ', '', seq[1])
        seqname = seq[0]
        seqname = re.sub(r'^>', '', seqname)
        name_short = re.search('^(\S+)', seqname).group(1)
        seqs = seq[1]
        seqdict[name_short] = seqs
    f.close()
    return seqdict


def get_mRNA_ranges(minimap_paf):
    """
    Get the mRNA ranges from a minimap2 PAF file.
    Parameters:
    - minimap_paf: Path to the minimap2 PAF file
    Returns:
    - A list of tuples containing the chromosome, start, and end positions of each mRNA range
    """
    mRNA_ranges = []
    with open(minimap_paf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[5]
            start = int(fields[7])
            end = int(fields[8])
            strand = fields[4]
            mRNA_ranges.append((chrom, start, end, strand))
    return mRNA_ranges


def get_exonerate_gene_range(exonerate_gff):
    gene_range = []
    with open(exonerate_gff, 'r') as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) > 4 and fields[2] == 'gene':
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                gene_range.append((chrom, start-21, end+21, strand))
    return gene_range


def get_intron_ranges(exonerate_gff):
    """
    Get the intron ranges from an Exonerate GFF file.
    Parameters:
    - exonerate_gff: Path to the Exonerate GFF file
    Returns:
    - A list of tuples containing the chromosome, start, and end positions of each intron range
    """
    intron_ranges = []
    with open(exonerate_gff, 'r') as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) > 4 and fields[2] == 'intron':
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                intron_ranges.append((chrom, start, end, strand))
    return intron_ranges


def get_exon_range_from_mRNA_and_intron_ranges(mRNA_ranges, intron_ranges):
    mRNA_dict = {}
    num = 0
    for chrom, start, end, strand in mRNA_ranges:
        mRNA_id = f'mRNA_{num}'
        mRNA_dict[mRNA_id] = {'chrom': chrom, 'start': start,
                              'end': end, 'strand': strand, 'exons': [], 'introns': []}
        for intron_chrom, intron_start, intron_end, intron_strand in intron_ranges:
            if chrom == intron_chrom and strand == intron_strand:
                if section((start, end), (intron_start, intron_end), int_flag=True, just_judgement=True):
                    mRNA_dict[mRNA_id]['introns'].append(
                        (intron_start, intron_end))
        mRNA_dict[mRNA_id]['exons'] = interval_minus_set(
            (start, end), mRNA_dict[mRNA_id]['introns'])
        mRNA_dict[mRNA_id]['exons'] = sorted(
            mRNA_dict[mRNA_id]['exons'], key=lambda x: x[0])
        mRNA_dict[mRNA_id]['introns'] = sorted(
            mRNA_dict[mRNA_id]['introns'], key=lambda x: x[0])

    return mRNA_dict


def get_range_haplotype(chr_id, start, end, bam_file, genome_file, output_dir, debug=False):
    """
    Get the haplotype sequences of a specific region.
    Parameters:
    - chr_id: Chromosome name
    - start: Start position of the region
    - end: End position of the region
    - bam_file: Path to the original BAM file
    - genome_file: Path to the reference genome file
    - work_dir: Path to the working directory
    """
    mkdir(output_dir)
    output_assem_h1_file = f"{output_dir}/range_hap1.fasta"
    output_assem_h2_file = f"{output_dir}/range_hap2.fasta"
    output_ref_file = f"{output_dir}/range.ref.fa"

    if os.path.exists(output_assem_h1_file) and os.path.getsize(output_assem_h1_file) > 0:
        print(
            f"Output files already exist: {output_assem_h1_file}, {output_assem_h2_file}, skipping reassembly.")
        return output_ref_file, output_assem_h1_file, output_assem_h2_file

    tmp_dir = f"{output_dir}/tmp_range_haplotype"
    mkdir(tmp_dir)

    cmd_string = "samtools view -bS %s %s:%d-%d > range.bam" % (
        bam_file, chr_id, start, end)
    cmd_run(cmd_string, cwd=tmp_dir)

    cmd_string = "samtools index range.bam"
    cmd_run(cmd_string, cwd=tmp_dir)

    cmd_string = "freebayes -f %s range.bam > range_variants.vcf" % (
        genome_file)
    cmd_run(cmd_string, cwd=tmp_dir)

    cmd_string = "whatshap phase -o range_phased.vcf --reference=%s range_variants.vcf range.bam" % (
        genome_file)
    cmd_run(cmd_string, cwd=tmp_dir)

    cmd_string = "bgzip range_phased.vcf && tabix range_phased.vcf.gz"
    cmd_run(cmd_string, cwd=tmp_dir)

    cmd_string = "samtools faidx %s %s:%d-%d > range.ref.fa" % (
        genome_file, chr_id, start, end)
    cmd_run(cmd_string, cwd=tmp_dir)

    cmd_string = "bcftools consensus -H 1 -f range.ref.fa range_phased.vcf.gz > range_hap1.fasta" % (
    )
    cmd_run(cmd_string, cwd=tmp_dir)

    cmd_string = "bcftools consensus -H 2 -f range.ref.fa range_phased.vcf.gz > range_hap2.fasta" % (
    )
    cmd_run(cmd_string, cwd=tmp_dir)

    hap1_file = "%s/range_hap1.fasta" % (tmp_dir)
    hap2_file = "%s/range_hap2.fasta" % (tmp_dir)
    ref_file = "%s/range.ref.fa" % (tmp_dir)

    cmd_run(f"mv {hap1_file} {output_assem_h1_file}", cwd=tmp_dir)
    cmd_run(f"mv {hap2_file} {output_assem_h2_file}", cwd=tmp_dir)
    cmd_run(f"mv {ref_file} {output_ref_file}", cwd=tmp_dir)

    if debug is False:
        rmdir(tmp_dir)

    return output_ref_file, output_assem_h1_file, output_assem_h2_file


def get_range_assembly(chr_id, start, end, bam_file, genome_file, output_dir=None, debug=False, assembly_tool='spades', polish=True):
    """
    Get the assembly sequences of a specific region.
    Parameters:
    - chr_id: Chromosome name
    - start: Start position of the region
    - end: End position of the region
    - bam_file: Path to the original BAM file
    - genome_file: Path to the reference genome file
    - work_dir: Path to the working directory
    """
    mkdir(output_dir)
    output_ref_file = f"{output_dir}/range.ref.fa"
    output_assem_file = f"{output_dir}/range.assem.fa"
    tmp_dir = f"{output_dir}/tmp_range_assem"
    mkdir(tmp_dir, keep=False)

    if os.path.exists(output_ref_file) and os.path.getsize(output_ref_file) > 0 and os.path.exists(output_assem_file) and os.path.getsize(output_assem_file) > 0:
        print(
            f"Output files already exist: {output_assem_file}, skipping reassembly.")

    # 1. 提取高质量成对 reads（MQ ≥ 30，proper pair），输出为两个 fastq
    cmd_string = f"samtools view -u -f 3 -q 30 {bam_file} {chr_id}:{start}-{end} | samtools collate -Ou - | samtools fastq -1 read_1.fq -2 read_2.fq -0 /dev/null -s /dev/null -n - > /dev/null"
    cmd_run(cmd_string, cwd=tmp_dir)

    # 2. 提取该区域的参考序列（用于 trusted contig）
    cmd_string = f"samtools faidx {genome_file} {chr_id}:{start}-{end} > {output_ref_file}"
    cmd_run(cmd_string, cwd=tmp_dir)

    # 1.1 如果没有成对 reads，则直接返回
    if not os.path.exists(f"{tmp_dir}/read_1.fq") or not os.path.exists(f"{tmp_dir}/read_2.fq"):
        cmd_run(f"touch {output_assem_file}", cwd=tmp_dir)
        if debug is False:
            rmdir(tmp_dir)
        return output_ref_file, output_assem_file

    if os.path.getsize(f"{tmp_dir}/read_1.fq") == 0 or os.path.getsize(f"{tmp_dir}/read_2.fq") == 0:
        cmd_run(f"touch {output_assem_file}", cwd=tmp_dir)
        if debug is False:
            rmdir(tmp_dir)
        return output_ref_file, output_assem_file

    try:
        if assembly_tool == 'spades':
            # 3. 使用 SPAdes 进行引导式拼接（输入左右端 reads）
            cmd_string = "spades.py -t 1 -m 5 -1 read_1.fq -2 read_2.fq -o range_spades_out"
            cmd_run(cmd_string, cwd=tmp_dir)
            # cmd_string = "spades.py --only-assembler -1 read_1.fq -2 read_2.fq --trusted-contigs range.ref.fa -o range_spades_out"
            # cmd_run(cmd_string, cwd=work_dir)
            assem_outfasta = f"{tmp_dir}/range_spades_out/contigs.fasta"

        elif assembly_tool == 'megahit':
            cmd_string = "megahit -t 1 -m 5e+9 -1 read_1.fq -2 read_2.fq -o range_megahit_out"
            cmd_run(cmd_string, cwd=tmp_dir)
            assem_outfasta = f"{tmp_dir}/range_megahit_out/final.contigs.fa"

        if polish:
            # 4. 将 reads 回帖到 SPAdes 拼出来的 contig 上
            cmd_string = f"bwa index {assem_outfasta}"
            cmd_run(cmd_string, cwd=tmp_dir)
            cmd_string = f"bwa mem {assem_outfasta} read_1.fq read_2.fq | samtools sort -o aln.bam"
            cmd_run(cmd_string, cwd=tmp_dir)
            cmd_string = "samtools index aln.bam"
            cmd_run(cmd_string, cwd=tmp_dir)

            # 5. 使用 Pilon 进行拼接纠错
            cmd_string = f"pilon --genome {assem_outfasta} --frags aln.bam --output polished --outdir polished_dir --vcf"
            cmd_run(cmd_string, cwd=tmp_dir)
            assem_outfasta = f"{tmp_dir}/polished_dir/polished.fasta"

    except:
        assem_outfasta = f"{tmp_dir}/non_exist.fasta"

    if not os.path.exists(assem_outfasta) or os.path.getsize(assem_outfasta) == 0:
        cmd_run(f"touch {output_assem_file}", cwd=tmp_dir)
    else:
        seq_dict = read_fasta(assem_outfasta)

        with open(output_assem_file, 'w') as out_f:
            num = 0
            for seq_id in sorted(seq_dict.keys(), key=lambda x: len(seq_dict[x]), reverse=True):
                new_seq_id = f"contig_{num}"
                out_f.write(f">{new_seq_id}\n{seq_dict[seq_id]}\n")
                num += 1

    if debug is False:
        rmdir(tmp_dir)


def get_range_annotation(local_assem_fasta, ref_pt_fasta, ref_cDNA_fasta, results_json_file, work_dir, debug=False, blast_tool='blast'):
    """
    Get the annotation of a specific region.
    Parameters:
    - local_assem_fasta: Path to the local assembly FASTA file
    - ref_pt_fasta: Path to the reference protein FASTA file
    - ref_cDNA_fasta: Path to the reference cDNA FASTA file
    - work_dir: Path to the working directory
    """
    mkdir(work_dir)

    cmd_string = f"exonerate --model protein2genome --showtargetgff yes --showquerygff no --showalignment no --minintron 20 --percent 30 --bestn 1 --maxintron 20000 {ref_pt_fasta} {local_assem_fasta} > exonerate.gff"
    cmd_run(cmd_string, cwd=work_dir)
    cmd_string = f"minimap2 -x splice -uf --secondary=no {local_assem_fasta} {ref_cDNA_fasta} > exon.aln.paf"
    # print(cmd_string)
    cmd_run(cmd_string, cwd=work_dir)

    contigs_seq_dict = read_fasta(local_assem_fasta)
    ref_pt_seq_dict = read_fasta(ref_pt_fasta)
    ref_pt_seq_length = len(ref_pt_seq_dict[list(ref_pt_seq_dict.keys())[0]])

    minimap_paf = f"{work_dir}/exon.aln.paf"
    exonerate_gff = f"{work_dir}/exonerate.gff"

    mRNA_ranges = get_mRNA_ranges(minimap_paf)
    intron_ranges = get_intron_ranges(exonerate_gff)
    exonerate_gene_range = get_exonerate_gene_range(exonerate_gff)
    mRNA_ranges = mRNA_ranges + exonerate_gene_range
    mRNA_dict = get_exon_range_from_mRNA_and_intron_ranges(
        mRNA_ranges, intron_ranges)

    for mRNA_id in mRNA_dict.keys():
        chrom = mRNA_dict[mRNA_id]['chrom']
        strand = mRNA_dict[mRNA_id]['strand']
        exons = mRNA_dict[mRNA_id]['exons']
        exon_seq = ''
        for exon_start, exon_end in exons:
            exon_seq += contigs_seq_dict[chrom][exon_start:exon_end + 1]
        exon_seq_len = len(exon_seq)
        mRNA_dict[mRNA_id]['exon_seq'] = exon_seq

        # if strand == '-':
        #     exon_seq = exon_seq[::-1].translate(str.maketrans('ATCG', 'TAGC'))

        tmp_dir = os.path.join(work_dir, 'tmp_' + str(mRNA_id))
        mkdir(tmp_dir, keep=False)
        exon_seq_file = os.path.join(tmp_dir, f"{mRNA_id}.exon.fna")
        with open(exon_seq_file, 'w') as f:
            f.write(f">{mRNA_id}\n{exon_seq}\n")

        cmd_string = f"TransDecoder.LongOrfs -S -m 10 -t {exon_seq_file}"
        cmd_run(cmd_string, cwd=tmp_dir)

        if blast_tool == 'diamond':
            bls_file = f"{tmp_dir}/{mRNA_id}.diamond.blastp.out"
            cmd_string = f"diamond blastp --query {exon_seq_file}.transdecoder_dir/longest_orfs.pep --db {ref_pt_fasta} --outfmt 6 --max-target-seqs 1 --evalue 1e-5 > {bls_file}"
            cmd_run(cmd_string, cwd=tmp_dir)
        elif blast_tool == 'blast':
            bls_file = f"{tmp_dir}/{mRNA_id}.blastp.out"
            if not os.path.exists(f'{ref_pt_fasta}.pdb'):
                cmd_string = f"makeblastdb -in {ref_pt_fasta} -dbtype prot"
                cmd_run(cmd_string, cwd=tmp_dir)
            cmd_string = f"blastp -query {exon_seq_file}.transdecoder_dir/longest_orfs.pep -db {ref_pt_fasta} -outfmt 6 -max_target_seqs 1 -evalue 1e-5 > {bls_file}"
            cmd_run(cmd_string, cwd=tmp_dir)

        hits = {}
        with open(bls_file, 'r') as blastp_file:
            for line in blastp_file:
                fields = line.strip().split('\t')
                query_id = fields[0]
                evalue = float(fields[10])
                identity = float(fields[2])/100
                subject_start = int(fields[8])
                subject_end = int(fields[9])
                subject_aln_length = abs(subject_end - subject_start) + 1
                subject_coverage = subject_aln_length / ref_pt_seq_length
                hits[query_id] = {
                    'evalue': evalue,
                    'identity': identity,
                    'coverage': subject_coverage
                }

        if len(hits) == 0:
            continue

        best_hit = sorted(hits.items(), key=lambda x: (
            x[1]['evalue'], -x[1]['identity'], -x[1]['coverage']))[0]
        best_hit_id = best_hit[0]
        cds_dict = read_fasta(
            f"{tmp_dir}/{mRNA_id}.exon.fna.transdecoder_dir/longest_orfs.cds")
        cds_seq = cds_dict[best_hit[0]]
        pt_dict = read_fasta(
            f"{tmp_dir}/{mRNA_id}.exon.fna.transdecoder_dir/longest_orfs.pep")
        pt_seq = pt_dict[best_hit[0]]
        mRNA_dict[mRNA_id]['cds_seq'] = cds_seq
        mRNA_dict[mRNA_id]['pt_seq'] = pt_seq
        mRNA_dict[mRNA_id]['evalue'] = best_hit[1]['evalue']
        mRNA_dict[mRNA_id]['identity'] = best_hit[1]['identity']
        mRNA_dict[mRNA_id]['coverage'] = best_hit[1]['coverage']

        utr5_len = 0
        utr3_len = 0
        with open(f"{tmp_dir}/{mRNA_id}.exon.fna.transdecoder_dir/longest_orfs.cds", 'r') as cds_file:
            for line in cds_file:
                if line.startswith('>'):
                    if best_hit_id == line.strip().split()[0].replace('>', ''):
                        match = re.search("(\d+)-(\d+)\(([+-])\)", line)
                        if match:
                            cds_start = int(match.group(1))
                            cds_end = int(match.group(2))
                            cds_strand = match.group(3).strip()
                            if cds_strand == '+':
                                utr5_len = cds_start - 1
                                utr3_len = exon_seq_len - cds_end
                            else:
                                cds_start, cds_end = cds_end, cds_start
                                utr5_len = exon_seq_len - cds_end
                                utr3_len = cds_start - 1

        utr5_list = []
        rest_utr5_len = utr5_len
        if mRNA_dict[mRNA_id]['strand'] == '+':
            for exon_start, exon_end in mRNA_dict[mRNA_id]['exons']:
                if rest_utr5_len <= 0:
                    break
                exon_length = exon_end - exon_start + 1
                if rest_utr5_len >= exon_length:
                    utr5_list.append((exon_start, exon_end))
                    rest_utr5_len -= exon_length
                else:
                    utr5_list.append(
                        (exon_start, exon_start + rest_utr5_len - 1))
                    rest_utr5_len = 0
        else:
            for exon_start, exon_end in reversed(mRNA_dict[mRNA_id]['exons']):
                if rest_utr5_len <= 0:
                    break
                exon_length = exon_end - exon_start + 1
                if rest_utr5_len >= exon_length:
                    utr5_list.append((exon_start, exon_end))
                    rest_utr5_len -= exon_length
                else:
                    utr5_list.append(
                        (exon_start + exon_length - rest_utr5_len, exon_end))
                    rest_utr5_len = 0

        utr3_list = []
        rest_utr3_len = utr3_len
        if mRNA_dict[mRNA_id]['strand'] == '+':
            for exon_start, exon_end in reversed(mRNA_dict[mRNA_id]['exons']):
                if rest_utr3_len <= 0:
                    break
                exon_length = exon_end - exon_start + 1
                if rest_utr3_len >= exon_length:
                    utr3_list.append((exon_start, exon_end))
                    rest_utr3_len -= exon_length
                else:
                    utr3_list.append((exon_end - rest_utr3_len + 1, exon_end))
                    rest_utr3_len = 0
        else:
            for exon_start, exon_end in mRNA_dict[mRNA_id]['exons']:
                if rest_utr3_len <= 0:
                    break
                exon_length = exon_end - exon_start + 1
                if rest_utr3_len >= exon_length:
                    utr3_list.append((exon_start, exon_end))
                    rest_utr3_len -= exon_length
                else:
                    utr3_list.append(
                        (exon_start, exon_start + rest_utr3_len - 1))
                    rest_utr3_len = 0

        mRNA_dict[mRNA_id]['utr5'] = utr5_list
        mRNA_dict[mRNA_id]['utr5'] = sorted(
            mRNA_dict[mRNA_id]['utr5'], key=lambda x: x[0])
        mRNA_dict[mRNA_id]['utr3'] = utr3_list
        mRNA_dict[mRNA_id]['utr3'] = sorted(
            mRNA_dict[mRNA_id]['utr3'], key=lambda x: x[0])
        mRNA_dict[mRNA_id]['cds'] = interval_minus_set(
            (mRNA_dict[mRNA_id]['start'], mRNA_dict[mRNA_id]['end']), utr5_list + utr3_list + mRNA_dict[mRNA_id]['introns'])

        # check_cds_seq = ''
        # for cds_start, cds_end in mRNA_dict[mRNA_id]['cds']:
        #     check_cds_seq += contigs_seq_dict[chrom][cds_start:cds_end + 1]
        # if mRNA_dict[mRNA_id]['strand'] == '-':
        #     check_cds_seq = check_cds_seq[::-
        #                                   1].translate(str.maketrans('ATCG', 'TAGC'))
        # if check_cds_seq != cds_seq:
        #     raise ValueError(f"CDS sequence mismatch")

        if debug is False:
            rmdir(tmp_dir)

    results_dict = {}
    for mRNA_id, mRNA_info in mRNA_dict.items():
        chrom = mRNA_info['chrom']
        results_dict.setdefault(
            chrom, {"seq": contigs_seq_dict[chrom], "mRNAs": {}})
        results_dict[chrom]["mRNAs"][mRNA_id] = mRNA_info

    for chrom in contigs_seq_dict:
        if chrom not in results_dict:
            results_dict[chrom] = {"seq": contigs_seq_dict[chrom], "mRNAs": {}}

    with open(results_json_file, 'w') as f:
        json.dump(results_dict, f, indent=4)

    if debug is False:
        rmdir(work_dir)

    return results_json_file


def get_cds_length(mRNA):
    cds_length = 0
    for cds in mRNA.sub_features:
        if cds.type == 'CDS':
            cds_length += abs(cds.start - cds.end) + 1
    return cds_length


def get_model_mRNA(gene):
    mRNA_id_list = [i for i in gene.sub_features if i.type == 'mRNA']
    model_mRNA = sorted(
        mRNA_id_list, key=lambda x: get_cds_length(x), reverse=True)[0]
    return model_mRNA


# local reassembly for genic region
def build_gene_db(genome_file, gene_gff_file, db_path, gene_flank=2000, intron_flank=500):
    mkdir(db_path)

    gene_dict = read_gff_file(gene_gff_file)['gene']
    contigs_seq_dict = read_fasta(genome_file)

    for g_id in gene_dict:
        gene_dir = os.path.join(db_path, g_id)
        mkdir(gene_dir)

        g = gene_dict[g_id]
        m = get_model_mRNA(g)
        exon_list = []
        for exon in m.sub_features:
            if exon.type == 'exon':
                exon_list.append((exon.start, exon.end))
        start_exon = sorted(exon_list, key=lambda x: x[0])[0]
        end_exon = sorted(exon_list, key=lambda x: x[1], reverse=True)[0]

        # get range list
        range_list = []
        if start_exon == end_exon:
            range_list.append(
                (start_exon[0] - gene_flank, end_exon[1] + gene_flank))
        else:
            range_list.append(
                (start_exon[0] - gene_flank, start_exon[1] + intron_flank))
            range_list.append(
                (end_exon[0] - intron_flank, end_exon[1] + gene_flank))
            for i in exon_list:
                if i == start_exon or i == end_exon:
                    continue
                range_list.append((i[0] - intron_flank, i[1] + intron_flank))

        range_list = sorted(range_list, key=lambda x: x[0])
        range_list = merge_intervals(range_list, int=True)
        range_list = [(max(0, r[0]), r[1])
                      for r in range_list]  # Ensure no negative start positions

        range_dict = {
            'gene_id': g_id,
            'chr_id': g.chr_id,
            'strand': m.strand,
            'range_list': range_list,
            'exon_list': exon_list,
        }

        with open(os.path.join(gene_dir, f"{g_id}.range.json"), 'w') as f:
            json.dump(range_dict, f, indent=4)

        # extract gene seq
        gene_seq = contigs_seq_dict[g.chr_id][
            g.start - gene_flank - 1:g.end + gene_flank]
        if m.strand == '-':
            gene_seq = gene_seq[::-1].translate(str.maketrans('ATCG', 'TAGC'))
        with open(os.path.join(gene_dir, f"{g_id}.gene.fa"), 'w') as f:
            f.write(f">{g_id}\n{gene_seq}\n")

        # extract exon seq
        exon_seq = ''
        for exon_start, exon_end in exon_list:
            exon_seq += contigs_seq_dict[g.chr_id][exon_start - 1:exon_end]
        if m.strand == '-':
            exon_seq = exon_seq[::-1].translate(str.maketrans('ATCG', 'TAGC'))

        exon_seq_file = os.path.join(gene_dir, f"{g_id}.exon.fa")
        with open(exon_seq_file, 'w') as f:
            f.write(f">{g_id}\n{exon_seq}\n")

        # extract aa seq
        cds_list = [cds for cds in m.sub_features if cds.type == 'CDS']
        if m.strand == '+':
            cds_list = sorted(cds_list, key=lambda x: x.start)
        else:
            cds_list = sorted(cds_list, key=lambda x: x.end, reverse=True)

        cds_seq = ''
        for cds in cds_list:
            if m.strand == '+':
                cds_seq += contigs_seq_dict[g.chr_id][cds.start - 1:cds.end]
            else:
                cds_seq += contigs_seq_dict[g.chr_id][cds.start -
                                                      1:cds.end][::-1].translate(str.maketrans('ATCG', 'TAGC'))

        aa_seq = cds_judgment(cds_seq, parse_phase=False,
                              keep_stop=True, return_cds=False)[2]

        with open(os.path.join(gene_dir, f"{g_id}.aa.fa"), 'w') as f:
            f.write(f">{g_id}\n{aa_seq}\n")

        with open(os.path.join(gene_dir, f"{g_id}.cds.fa"), 'w') as f:
            f.write(f">{g_id}\n{cds_seq}\n")


def gene_pipeline(gene_id, genome_file, gene_db_path, bam_file, work_dir=None, debug=False, assembly_mode='assembly', assembly_tool='spades', polish=True, remove_existing=False):
    mkdir(work_dir, keep=not remove_existing)

    final_results_json_file = os.path.join(work_dir, f"anno.json")
    final_stitch_fasta = os.path.join(work_dir, f"stitch.exon.assem.fa")
    if os.path.exists(final_results_json_file) and os.path.getsize(final_results_json_file) > 0 and os.path.exists(final_stitch_fasta) and os.path.getsize(final_stitch_fasta) > 0:
        print(
            f"Output files already exist: {final_results_json_file}, {final_stitch_fasta}, skipping reassembly.")
        return final_results_json_file, final_stitch_fasta

    tmp_dir = os.path.join(work_dir, 'tmp_gene_pipeline')
    mkdir(tmp_dir, keep=False)

    assem_dir = os.path.join(tmp_dir, 'assem')
    mkdir(assem_dir)
    range_dict = json.load(
        open(os.path.join(gene_db_path, gene_id, f"{gene_id}.range.json"), 'r'))

    chr_id = range_dict['chr_id']
    strand = range_dict['strand']
    range_list = range_dict['range_list']
    # Sort by start position
    range_list = sorted(range_list, key=lambda x: x[0])
    range_dict = {f"range_{i}": (start, end)
                  for i, (start, end) in enumerate(range_list)}

    for range_key, (start, end) in range_dict.items():
        if assembly_mode == 'assembly':
            get_range_assembly(chr_id, start, end, bam_file, genome_file,
                               output_dir=f"{assem_dir}/{range_key}", debug=debug, assembly_tool=assembly_tool, polish=polish)
            assem_file_name = 'range.assem.fa'

    # Combine all assembly results
    exon_assem_fasta = os.path.join(assem_dir, 'range.exon.assem.fa')
    num = 0
    with open(exon_assem_fasta, 'w') as out_f:
        for range_key in sorted(range_dict.keys(), key=lambda x: range_dict[x][0]):
            range_assem_file = os.path.join(
                assem_dir, range_key, assem_file_name)
            range_assem_seq_dict = read_fasta(range_assem_file)
            for seq_id, seq in range_assem_seq_dict.items():
                new_seq_id = f"contig_{num}"
                out_f.write(f">{new_seq_id}\n{seq}\n")
                num += 1

    # blastn to ref gene seq
    ref_gene_seq_file = os.path.join(
        gene_db_path, gene_id, f"{gene_id}.gene.fa")
    cmd_run(f"cp {ref_gene_seq_file} {assem_dir}/ref.gene.fa", cwd=assem_dir)
    cmd_string = f"blastn -query {exon_assem_fasta} -subject ref.gene.fa -outfmt 6 -max_hsps 1 -max_target_seqs 1 -evalue 1e-5 > exon.blastn.out"
    cmd_run(cmd_string, cwd=assem_dir)

    blast_result_file = os.path.join(assem_dir, 'exon.blastn.out')
    blast_results = {}
    with open(blast_result_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            query_id = fields[0]
            sstart = int(fields[8])
            send = int(fields[9])
            strand = '+' if sstart < send else '-'
            if strand == '-':
                sstart, send = send, sstart
            blast_results[query_id] = (sstart, send, strand)

    exon_assem_dict = read_fasta(exon_assem_fasta)
    stitch_seq = ''
    assem_id_list = blast_results.keys()
    assem_id_list = sorted(assem_id_list, key=lambda x: blast_results[x][0])
    for assem_id in assem_id_list:
        sstart, send, strand = blast_results[assem_id]
        seq = exon_assem_dict[assem_id]
        if strand == '-':
            seq = seq[::-1].translate(str.maketrans('ATCG', 'TAGC'))
        stitch_seq += seq + 'N' * 100

    stitch_fasta = os.path.join(assem_dir, 'stitch.exon.assem.fa')
    with open(stitch_fasta, 'w') as f:
        f.write(f">stitch_exon_assem\n{stitch_seq}\n")

    # run annotation
    ref_pt_fasta = os.path.join(gene_db_path, gene_id, f"{gene_id}.aa.fa")
    ref_cDNA_fasta = os.path.join(gene_db_path, gene_id, f"{gene_id}.exon.fa")
    anno_dir = os.path.join(tmp_dir, 'anno')
    mkdir(anno_dir)
    results_json_file = os.path.join(anno_dir, f"anno.json")
    get_range_annotation(stitch_fasta, ref_pt_fasta, ref_cDNA_fasta,
                         results_json_file, anno_dir+"/tmp", debug=debug)

    cmd_run(f"cp {results_json_file} {work_dir}/", cwd=anno_dir)
    cmd_run(f"cp {stitch_fasta} {work_dir}/", cwd=anno_dir)

    if debug is False:
        rmdir(tmp_dir)
        shutil.make_archive(work_dir, 'zip', work_dir)
        rmdir(work_dir)


def time_now():
    time_tmp = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    return time_tmp


def abortable_worker(func, *args, **kwargs):
    timeout = kwargs.get('timeout', None)
    p = ThreadPool(1)
    res = p.apply_async(func, args=args)
    try:
        out = res.get(timeout)  # Wait timeout seconds for func to complete.
        return out
    except TimeoutError:
        return "Aborting due to timeout"
        raise


def logging_init(program_name, log_file=None, log_level=logging.DEBUG, console_level=logging.ERROR):
    # create logger with 'program_name'
    logger = logging.getLogger(program_name)
    logger.setLevel(logging.DEBUG)
    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    if not log_file is None:
        # create file handler which logs even debug messages
        fh = logging.FileHandler(log_file)
        fh.setLevel(log_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(console_level)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(ch)

    return logger


def get_more_para(f, para_tuple):
    return f(*para_tuple)


def multiprocess_running(f, args_dict, pool_num, log_file=None, silence=False, timeout=None):
    """
    :param f: function name
    :param args_dict: a dict, value is every args in a tuple, key is the job_id
    :param pool_num: process num of running at same time
    :param silence:
    :return:

    example

    def f(x, y):
        time.sleep(1)
        return x * y

    args_list = list(zip(range(1,200),range(2,201)))
    pool_num = 5

    multiprocess_running(f, args_list, pool_num, log_file='/lustre/home/xuyuxing/Work/Other/saif/Meth/tmp/log')

    """
    num_tasks = len(list(args_dict))

    module_log = logging_init(f.__name__, log_file)
    if not silence:
        print('args_list have %d object and pool_num is %d' %
              (num_tasks, pool_num))
    module_log.info('running with mulitprocess')
    module_log.info('args_list have %d object and pool_num is %d' %
                    (num_tasks, pool_num))

    p_dict = {}

    args_id_list = list(args_dict.keys())
    args_list = [args_dict[i] for i in args_id_list]

    if len(args_dict) > 0:
        f_with_para = partial(get_more_para, f)

        start_time = time.time()
        if not silence:
            print(time_now() + '\tBegin: ')
        module_log.info('Begin: ')

        if timeout is None:
            with Pool(processes=pool_num) as pool:
                # for i, output in enumerate(pool.imap_unordered(f_with_para, args_list, chunksize=1)):
                for i, output in enumerate(pool.imap(f_with_para, args_list, chunksize=1)):
                    job_id = args_id_list[i]

                    p_dict[job_id] = {
                        'args': args_list[i],
                        'output': output,
                        'error': None
                    }
                    # print(i)
                    round_time = time.time()
                    if round_time - start_time > 5:
                        if not silence:
                            print(time_now() + '\t%d/%d %.2f%% parsed' %
                                  (i, num_tasks, i / num_tasks * 100))
                        module_log.info('%d/%d %.2f%% parsed' %
                                        (i, num_tasks, i / num_tasks * 100))
                        start_time = round_time
                    # module_log.info('%d/%d %.2f%% parsed' % (i, num_tasks, i / num_tasks * 100))
        else:
            abortable_func = partial(
                abortable_worker, f_with_para, timeout=timeout)

            # with Pool(processes=pool_num, maxtasksperchild=1) as pool:
            with Pool(processes=pool_num) as pool:
                # for i, output in enumerate(pool.imap_unordered(f_with_para, args_list, chunksize=1)):
                it = pool.imap(abortable_func, args_list, chunksize=1)
                i = -1
                while 1:
                    i += 1
                    try:
                        output = it.next()
                        if args_id_list:
                            job_id = args_id_list[i]
                        else:
                            job_id = 'ID_' + str(i)

                        if output == "Aborting due to timeout":
                            p_dict[job_id] = {
                                'args': args_list[i],
                                'output': None,
                                'error': "timeout"
                            }
                        else:
                            p_dict[job_id] = {
                                'args': args_list[i],
                                'output': output,
                                'error': None
                            }

                        round_time = time.time()
                        if round_time - start_time > 5:
                            if not silence:
                                print(time_now() + '\t%d/%d %.2f%% parsed' %
                                      (i, num_tasks, i / num_tasks * 100))
                            module_log.info('%d/%d %.2f%% parsed' %
                                            (i, num_tasks, i / num_tasks * 100))
                            start_time = round_time

                    except StopIteration:
                        break

        module_log.info('%d/%d %.2f%% parsed' %
                        (i, num_tasks, i / num_tasks * 100))
        module_log.info('All args_list task finished')

        if not silence:
            print(time_now() + '\t%d/%d %.2f%% parsed' %
                  (i, num_tasks, i / num_tasks * 100))
            print(time_now() + '\tAll args_list task finished')

    del module_log.handlers[:]

    return p_dict


def genome_pipeline(genome_file, gene_db_path, bam_file, work_dir=None, debug=False, assembly_mode='assembly', assembly_tool='spades', polish=True, threads=20):
    mkdir(work_dir)
    gene_id_list = os.listdir(gene_db_path)

    ok_gene_id_list = os.listdir(work_dir)
    ok_gene_id_list = [i for i in ok_gene_id_list if i.endswith('.zip')]
    ok_gene_id_list = [i.replace('.zip', '') for i in ok_gene_id_list]

    gene_id_list = list(set(gene_id_list) - set(ok_gene_id_list))

    args_dict = {}
    for gene_id in gene_id_list:
        gene_dir = os.path.join(gene_db_path, gene_id)
        if not os.path.exists(os.path.join(gene_dir, f"{gene_id}.range.json")):
            continue
        args_dict[gene_id] = (gene_id, genome_file, gene_db_path,
                              bam_file, os.path.join(work_dir, gene_id), debug, assembly_mode, assembly_tool, polish, True)

    mlt_dict = multiprocess_running(
        gene_pipeline, args_dict, threads, log_file=os.path.join(work_dir, "gene_pipeline.log"))
