import subprocess
import os
import shutil
import re
import json


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


def get_range_haplotype(chr_id, start, end, bam_file, genome_file, output_ref_file, output_assem_h1_file, output_assem_h2_file, work_dir, debug=False, return_ref=False):
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

    if os.path.exists(output_assem_h1_file) and os.path.getsize(output_assem_h1_file) > 0:
        print(
            f"Output files already exist: {output_assem_h1_file}, {output_assem_h2_file}, skipping reassembly.")
        if return_ref:
            return output_ref_file, output_assem_h1_file, output_assem_h2_file
        else:
            return output_assem_h1_file, output_assem_h2_file

    mkdir(work_dir)

    cmd_string = "samtools view -bS %s %s:%d-%d > range.bam" % (
        bam_file, chr_id, start, end)
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "samtools index range.bam"
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "freebayes -f %s range.bam > range_variants.vcf" % (
        genome_file)
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "whatshap phase -o range_phased.vcf --reference=%s range_variants.vcf range.bam" % (
        genome_file)
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "bgzip range_phased.vcf && tabix range_phased.vcf.gz"
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "samtools faidx %s %s:%d-%d > range.ref.fa" % (
        genome_file, chr_id, start, end)
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "bcftools consensus -H 1 -f range.ref.fa range_phased.vcf.gz > range_hap1.fasta" % (
    )
    cmd_run(cmd_string, cwd=work_dir)

    cmd_string = "bcftools consensus -H 2 -f range.ref.fa range_phased.vcf.gz > range_hap2.fasta" % (
    )
    cmd_run(cmd_string, cwd=work_dir)

    hap1_file = "%s/range_hap1.fasta" % (work_dir)
    hap2_file = "%s/range_hap2.fasta" % (work_dir)
    ref_file = "%s/range.ref.fa" % (work_dir)

    cmd_run(f"mv {hap1_file} {output_assem_h1_file}", cwd=work_dir)
    cmd_run(f"mv {hap2_file} {output_assem_h2_file}", cwd=work_dir)
    cmd_run(f"mv {ref_file} {output_ref_file}", cwd=work_dir)

    if debug is False:
        rmdir(work_dir)

    if return_ref:
        return output_ref_file, output_assem_h1_file, output_assem_h2_file
    else:
        return output_assem_h1_file, output_assem_h2_file


def get_range_assembly(chr_id, start, end, bam_file, genome_file, output_ref_file, output_assem_file, work_dir, debug=False, return_ref=False):
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
    if os.path.exists(output_ref_file) and os.path.getsize(output_ref_file) > 0 and os.path.exists(output_assem_file) and os.path.getsize(output_assem_file) > 0:
        print(
            f"Output files already exist: {output_assem_file}, skipping reassembly.")
        if return_ref:
            return output_ref_file, output_assem_file
        else:
            return output_assem_file

    mkdir(work_dir)

    # 1. 提取高质量成对 reads（MQ ≥ 30，proper pair），输出为两个 fastq
    cmd_string = f"samtools view -u -f 3 -q 30 {bam_file} {chr_id}:{start}-{end} | samtools collate -Ou - | samtools fastq -1 read_1.fq -2 read_2.fq -0 /dev/null -s /dev/null -n - > /dev/null"
    cmd_run(cmd_string, cwd=work_dir)

    # 2. 提取该区域的参考序列（用于 trusted contig）
    cmd_string = f"samtools faidx {genome_file} {chr_id}:{start}-{end} > {output_ref_file}"
    cmd_run(cmd_string, cwd=work_dir)

    # 1.1 如果没有成对 reads，则直接返回
    if not os.path.exists(f"{work_dir}/read_1.fq") or not os.path.exists(f"{work_dir}/read_2.fq"):
        cmd_run(f"touch {output_assem_file}", cwd=work_dir)
        if debug is False:
            rmdir(work_dir)
        return output_ref_file, output_assem_file

    if os.path.getsize(f"{work_dir}/read_1.fq") == 0 or os.path.getsize(f"{work_dir}/read_2.fq") == 0:
        cmd_run(f"touch {output_assem_file}", cwd=work_dir)
        if debug is False:
            rmdir(work_dir)
        return output_ref_file, output_assem_file

    # 3. 使用 SPAdes 进行引导式拼接（输入左右端 reads）
    cmd_string = "spades.py -1 read_1.fq -2 read_2.fq -o range_spades_out"
    cmd_run(cmd_string, cwd=work_dir)

    # cmd_string = "spades.py --only-assembler -1 read_1.fq -2 read_2.fq --trusted-contigs range.ref.fa -o range_spades_out"
    # cmd_run(cmd_string, cwd=work_dir)

    # 4. 将 reads 回帖到 SPAdes 拼出来的 contig 上
    cmd_string = "bwa index range_spades_out/contigs.fasta"
    cmd_run(cmd_string, cwd=work_dir)
    cmd_string = "bwa mem range_spades_out/contigs.fasta read_1.fq read_2.fq | samtools sort -o aln.bam"
    cmd_run(cmd_string, cwd=work_dir)
    cmd_string = "samtools index aln.bam"
    cmd_run(cmd_string, cwd=work_dir)

    # 5. 使用 Pilon 进行拼接纠错
    cmd_string = "pilon --genome range_spades_out/contigs.fasta --frags aln.bam --output polished --outdir polished_dir --vcf"
    cmd_run(cmd_string, cwd=work_dir)

    # 6. 将纠错后的 contig 提取出来
    polished_fasta = f"{work_dir}/polished_dir/polished.fasta"
    polished_seq_dict = read_fasta(polished_fasta)

    with open(output_assem_file, 'w') as out_f:
        num = 0
        for seq_id in sorted(polished_seq_dict.keys(), key=lambda x: len(polished_seq_dict[x]), reverse=True):
            new_seq_id = f"contig_{num}"
            out_f.write(f">{new_seq_id}\n{polished_seq_dict[seq_id]}\n")
            num += 1

    if debug is False:
        rmdir(work_dir)

    if return_ref:
        return output_ref_file, output_assem_file
    else:
        return output_assem_file


def get_range_annotation(local_assem_fasta, ref_pt_fasta, ref_cDNA_fasta, results_json_file, work_dir, debug=False):
    """
    Get the annotation of a specific region.
    Parameters:
    - local_assem_fasta: Path to the local assembly FASTA file
    - ref_pt_fasta: Path to the reference protein FASTA file
    - ref_cDNA_fasta: Path to the reference cDNA FASTA file
    - work_dir: Path to the working directory
    """
    if os.path.exists(results_json_file) and os.path.getsize(results_json_file) > 0:
        print(
            f"Results already exist: {results_json_file}, skipping reannotation.")
        return results_json_file

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

        cmd_string = f"TransDecoder.LongOrfs -t {exon_seq_file}"
        cmd_run(cmd_string, cwd=tmp_dir)

        diamond_bls_file = f"{tmp_dir}/{mRNA_id}.blastp.out"
        cmd_string = f"diamond blastp --query {exon_seq_file}.transdecoder_dir/longest_orfs.pep --db {ref_pt_fasta} --outfmt 6 --max-target-seqs 1 --evalue 1e-5 > {diamond_bls_file}"
        cmd_run(cmd_string, cwd=tmp_dir)

        hits = {}
        with open(diamond_bls_file, 'r') as blastp_file:
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

        check_cds_seq = ''
        for cds_start, cds_end in mRNA_dict[mRNA_id]['cds']:
            check_cds_seq += contigs_seq_dict[chrom][cds_start:cds_end + 1]
        if mRNA_dict[mRNA_id]['strand'] == '-':
            check_cds_seq = check_cds_seq[::-
                                          1].translate(str.maketrans('ATCG', 'TAGC'))
        if check_cds_seq != cds_seq:
            raise ValueError(f"CDS sequence mismatch")

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

    if debug is False:
        rmdir(work_dir)

    with open(results_json_file, 'w') as f:
        json.dump(results_dict, f, indent=4)

    return results_json_file