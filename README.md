# Local Reassembly

Local Reassembly is a tool for locally reassembling reads from a mapped BAM file. We have two modes of operation: "assembly" and "haplotype". In "assembly" mode, we use SPAdes to perform de novo assembly of the mapped reads. In "haplotype" mode, we use whatshap to perform haplotype assembly of the mapped reads. The output is a FASTA file containing the assembled contigs.

## Installation

Local Reassembly is a Python package that can be installed using pip. To install the package, run the following command:

```bash
pip install local_reassembly
```

We also recommend installing with Docker, please see Dockerfile for more information.

## Usage

- Assembly

```bash
usage: locasm reassm [-h] [-o OUTPUT_PREFIX] [-t TMP_WORK_DIR] [-d] [-m {assembly,haplotype}] input_genome_file input_bam_file region

Local reassembly

positional arguments:
  input_genome_file     input genome file in FASTA format
  input_bam_file        input BAM file
  region                genomic region in the format chr:start-end

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        output prefix for the reassembly files
  -t TMP_WORK_DIR, --tmp_work_dir TMP_WORK_DIR
                        temporary working directory, default is current directory
  -d, --debug           debug mode, default False
  -m {assembly,haplotype}, --mode {assembly,haplotype}
                        mode of operation: "assembly" for local assembly, "haplotype" for haplotype reconstruction
```

- Annotation

```bash
usage: locasm reanno [-h] [-o OUTPUT_PREFIX] [-t TMP_WORK_DIR] [-d] local_assem_fasta ref_pt_fasta ref_cDNA_fasta

Local reannotation

positional arguments:
  local_assem_fasta     local assembly FASTA file
  ref_pt_fasta          reference point FASTA file
  ref_cDNA_fasta        reference cDNA FASTA file

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        output prefix for the reannotation files
  -t TMP_WORK_DIR, --tmp_work_dir TMP_WORK_DIR
                        temporary working directory, default is current directory
  -d, --debug           debug mode, default False
```


