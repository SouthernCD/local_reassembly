import argparse
import uuid
import os
from local_reassembly.src import get_range_haplotype, get_range_assembly, get_range_annotation, build_gene_db, gene_pipeline


class CustomHelpFormatter(argparse.HelpFormatter):
    def add_subparsers(self, *args, **kwargs):
        subparsers_action = super().add_subparsers(*args, **kwargs)
        subparsers_action._parser_class = CustomSubcommandParser
        return subparsers_action


class CustomSubcommandParser(argparse.ArgumentParser):
    def format_help(self):
        formatter = self._get_formatter()

        # Add the usage
        formatter.add_usage(self.usage, self._actions,
                            self._mutually_exclusive_groups)

        # Add the description
        formatter.add_text(self.description)

        # Add the subcommands
        for action in self._actions:
            if isinstance(action, argparse._SubParsersAction):
                formatter.start_section("subcommands")
                for choice, subparser in action.choices.items():
                    formatter.add_text(f"{choice}: {subparser.description}\n")
                formatter.end_section()

        # Add the epilog
        formatter.add_text(self.epilog)

        # Return the full help string
        return formatter.format_help()


class Job(object):
    def __init__(self):
        pass

    def run_arg_parser(self):
        # argument parse

        parser = argparse.ArgumentParser(
            prog='reloc',
            description="Local reassembly and reannotation tool",
            formatter_class=CustomHelpFormatter
        )

        subparsers = parser.add_subparsers(
            title='subcommands', dest="subcommand_name")

        # argparse for reassm
        parser_a = subparsers.add_parser('reassm',
                                         description='Local reassembly',
                                         help='Local reassembly')
        parser_a.add_argument('input_genome_file', type=str,
                              help='input genome file in FASTA format')
        parser_a.add_argument('input_bam_file', type=str,
                              help='input BAM file')
        parser_a.add_argument('region', type=str,
                              help='genomic region in the format chr:start-end')
        parser_a.add_argument('-o', '--output_dir', type=str,
                              help='output directory for the reassembly files',
                              default='reassm_output')
        parser_a.add_argument('-d', '--debug', action='store_true',
                              help='debug mode, default False')
        parser_a.add_argument('-m', '--mode', type=str, choices=['assembly', 'haplotype'],
                              help='mode of operation: "assembly" for local assembly, "haplotype" for haplotype reconstruction',
                              default='assembly')
        parser_a.add_argument('-p', '--polish', action='store_true',
                              help='whether to polish the assembly with Pilon, default False')
        parser_a.add_argument('-a', '--assembler', type=str, choices=['spades', 'megahit'],
                              help='assembler to use, default is megahit',
                              default='megahit')

        # argparse for reanno
        parser_b = subparsers.add_parser('reanno',
                                         description='Local reannotation',
                                         help='Local reannotation')
        parser_b.add_argument('local_assem_fasta', type=str,
                              help='local assembly FASTA file')
        parser_b.add_argument('ref_pt_fasta', type=str,
                              help='reference point FASTA file')
        parser_b.add_argument('ref_cDNA_fasta', type=str,
                              help='reference cDNA FASTA file')
        parser_b.add_argument('-o', '--output_prefix', type=str,
                              help='output prefix for the reannotation files',
                              default='reanno_output')
        parser_b.add_argument('-t', '--tmp_work_dir', type=str,
                              help='temporary working directory, default is current directory',
                              default=None)
        parser_b.add_argument('-d', '--debug', action='store_true',
                              help='debug mode, default False')

        # argparse for genedb
        parser_c = subparsers.add_parser('genedb',
                                         description='Gene database generation',
                                         help='Gene database generation')
        parser_c.add_argument('input_genome_file', type=str,
                              help='input genome file in FASTA format')
        parser_c.add_argument('gene_gff_file', type=str,
                              help='gene annotation GFF file')
        parser_c.add_argument('db_path', type=str,
                              help='output database path')
        parser_c.add_argument('-g', '--gene_flank', type=int,
                              help='gene flanking region size, default 2000',
                              default=2000)
        parser_c.add_argument('-i', '--intron_flank', type=int,
                              help='intron flanking region size, default 500',
                              default=500)

        # argparse for gene pipeline
        parser_d = subparsers.add_parser('genepipe',
                                         description='Gene pipeline',
                                         help='Gene pipeline')
        parser_d.add_argument('gene_id', type=str,
                              help='gene ID to process')
        parser_d.add_argument('genome_file', type=str,
                              help='input genome file in FASTA format')
        parser_d.add_argument('gene_db_path', type=str,
                              help='path to the gene database, should be generated by genedb command')
        parser_d.add_argument('bam_file', type=str,
                              help='input BAM file')
        parser_d.add_argument('-w', '--work_dir', type=str,
                              help='working directory, default is current directory',
                              default=None)
        parser_d.add_argument('-d', '--debug', action='store_true',
                              help='debug mode, default False')
        parser_d.add_argument('-m', '--assembly_mode', type=str,
                              help='assembly mode, default is assembly',
                              default='assembly')
        parser_d.add_argument('-a', '--assembler', type=str, choices=['spades', 'megahit'],
                              help='assembler to use, default is megahit',
                              default='megahit')
        parser_d.add_argument('-p', '--polish', action='store_true',
                              help='whether to polish the assembly with Pilon, default False')

        self.arg_parser = parser

        self.args = parser.parse_args()

    def run(self):
        self.run_arg_parser()

        if self.args.subcommand_name == 'reassm':
            # Parse the region
            region = self.args.region
            if ':' not in region or '-' not in region:
                raise ValueError("Region must be in the format chr:start-end")
            chr_id, pos_range = region.split(':')
            start, end = map(int, pos_range.split('-'))

            # Prepare file paths
            bam_file = self.args.input_bam_file
            bam_file = os.path.abspath(bam_file)
            genome_file = self.args.input_genome_file
            genome_file = os.path.abspath(genome_file)
            output_dir = self.args.output_dir
            if output_dir is None:
                output_dir = f'./reassm_{uuid.uuid4().hex}'
            output_dir = os.path.abspath(output_dir)
            debug = self.args.debug
            assembly_tool = self.args.assembler
            polish = self.args.polish

            # Run the appropriate function based on mode
            if self.args.mode == 'assembly':
                get_range_assembly(chr_id, start, end, bam_file, genome_file, output_dir=output_dir,
                                   debug=debug, assembly_tool=assembly_tool, polish=polish)
            elif self.args.mode == 'haplotype':
                get_range_haplotype(chr_id, start, end, bam_file, genome_file, output_ref_file,
                                    output_assem_h1_file, output_assem_h2_file, work_dir, debug=debug, return_ref=debug)
            else:
                raise ValueError(
                    "Invalid mode. Choose 'assembly' or 'haplotype'.")

        elif self.args.subcommand_name == 'reanno':
            local_assem_fasta = self.args.local_assem_fasta
            local_assem_fasta = os.path.abspath(local_assem_fasta)
            ref_pt_fasta = self.args.ref_pt_fasta
            ref_pt_fasta = os.path.abspath(ref_pt_fasta)
            ref_cDNA_fasta = self.args.ref_cDNA_fasta
            ref_cDNA_fasta = os.path.abspath(ref_cDNA_fasta)
            output_prefix = self.args.output_prefix
            output_json_file = f"{output_prefix}_annotation.json"
            output_json_file = os.path.abspath(output_json_file)
            work_dir = self.args.tmp_work_dir
            if work_dir is None:
                work_dir = f'./reanno_{uuid.uuid4().hex}'
            work_dir = os.path.abspath(work_dir)
            debug = self.args.debug
            get_range_annotation(local_assem_fasta, ref_pt_fasta,
                                 ref_cDNA_fasta, output_json_file, work_dir, debug=debug)

        elif self.args.subcommand_name == 'genedb':
            input_genome_file = self.args.input_genome_file
            input_genome_file = os.path.abspath(input_genome_file)
            gene_gff_file = self.args.gene_gff_file
            gene_gff_file = os.path.abspath(gene_gff_file)
            db_path = self.args.db_path
            db_path = os.path.abspath(db_path)
            gene_flank = self.args.gene_flank
            intron_flank = self.args.intron_flank

            build_gene_db(input_genome_file, gene_gff_file, db_path,
                          gene_flank=gene_flank, intron_flank=intron_flank)

        elif self.args.subcommand_name == 'genepipe':
            gene_id = self.args.gene_id
            genome_file = self.args.genome_file
            genome_file = os.path.abspath(genome_file)
            gene_db_path = self.args.gene_db_path
            gene_db_path = os.path.abspath(gene_db_path)
            bam_file = self.args.bam_file
            bam_file = os.path.abspath(bam_file)
            work_dir = self.args.work_dir
            if work_dir is None:
                work_dir = f'./{gene_id}'
            work_dir = os.path.abspath(work_dir)
            debug = self.args.debug
            assembly_mode = self.args.assembly_mode
            assembler = self.args.assembler
            polish = self.args.polish
            gene_pipeline(gene_id, genome_file, gene_db_path,
                          bam_file, work_dir=work_dir, debug=debug, assembly_mode=assembly_mode, assembly_tool=assembler, polish=polish)

        else:
            self.arg_parser.print_help()


def main():
    job = Job()
    job.run()


if __name__ == '__main__':
    main()
