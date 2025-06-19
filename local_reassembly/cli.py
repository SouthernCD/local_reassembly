import argparse
import uuid
import os
from local_reassembly.src import get_range_haplotype, get_range_assembly, get_range_annotation


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
            prog='locasm',
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
        parser_a.add_argument('-o', '--output_prefix', type=str,
                              help='output prefix for the reassembly files',
                              default='reassm_output')
        parser_a.add_argument('-t', '--tmp_work_dir', type=str,
                              help='temporary working directory, default is current directory',
                              default=None)
        parser_a.add_argument('-d', '--debug', action='store_true',
                              help='debug mode, default False')
        parser_a.add_argument('-m', '--mode', type=str, choices=['assembly', 'haplotype'],
                              help='mode of operation: "assembly" for local assembly, "haplotype" for haplotype reconstruction',
                              default='assembly')

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
            output_prefix = self.args.output_prefix
            work_dir = self.args.tmp_work_dir
            if work_dir is None:
                work_dir = f'./reassm_{uuid.uuid4().hex}'
            work_dir = os.path.abspath(work_dir)
            debug = self.args.debug

            # Output files
            output_ref_file = f"{output_prefix}_ref.fasta"
            output_assem_file = f"{output_prefix}_assembly.fasta"
            output_assem_h1_file = f"{output_prefix}_haplotype_h1.fasta"
            output_assem_h2_file = f"{output_prefix}_haplotype_h2.fasta"
            output_ref_file = os.path.abspath(output_ref_file)
            output_assem_file = os.path.abspath(output_assem_file)
            output_assem_h1_file = os.path.abspath(output_assem_h1_file)
            output_assem_h2_file = os.path.abspath(output_assem_h2_file)

            # Run the appropriate function based on mode
            if self.args.mode == 'assembly':
                get_range_assembly(chr_id, start, end, bam_file, genome_file,
                                   output_ref_file, output_assem_file, work_dir, debug=debug, return_ref=debug)
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
        else:
            self.arg_parser.print_help()


def main():
    job = Job()
    job.run()


if __name__ == '__main__':
    main()
