import subprocess
import os
import shutil

def mkdir(dir_name, keep=True):
    if keep is False:
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
        os.makedirs(dir_name)
    else:
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

    return dir_name


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


def get_range_haplotype(chr_id, start, end, bam_file, genome_file, work_dir):
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

    return hap1_file, hap2_file, ref_file
