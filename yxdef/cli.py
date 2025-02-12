import argparse
from yxdef.src.pipeline import subcmd1_main, subcmd2_main

class CustomHelpFormatter(argparse.HelpFormatter):
    def add_subparsers(self, *args, **kwargs):
        subparsers_action = super().add_subparsers(*args, **kwargs)
        subparsers_action._parser_class = CustomSubcommandParser
        return subparsers_action

class CustomSubcommandParser(argparse.ArgumentParser):
    def format_help(self):
        formatter = self._get_formatter()
        
        # Add the usage
        formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups)
        
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
            prog='yxdef',
            description="Main command description.",
            formatter_class=CustomHelpFormatter
        )

        subparsers = parser.add_subparsers(
            title='subcommands', dest="subcommand_name")

        # argparse for subcmd1
        parser_a = subparsers.add_parser('subcmd1',
                                         description='Description for subcmd1',
                                         help='Description for subcmd1')
        parser_a.add_argument('input', type=str,
                              help='input file')
        parser_a.add_argument('output', type=str,
                              help='output file')
        parser_a.add_argument('-t', '--threads', type=int,
                              help='number of threads, default 1', default=1)
        parser_a.add_argument('-d', '--dry_run', action='store_true',
                              help='dry run, default False')
        parser_a.set_defaults(func=subcmd1_main)

        # argparse for subcmd2
        parser_b = subparsers.add_parser('subcmd2',
                                         description='Description for subcmd2',
                                         help='Description for subcmd2')
        parser_b.add_argument('input', type=str,
                              help='input file')
        parser_b.add_argument('output', type=str,
                              help='output file')
        parser_b.add_argument('-t', '--threads', type=int,
                              help='number of threads, default 1', default=1)
        parser_b.add_argument('-d', '--dry_run', action='store_true',
                              help='dry run, default False')
        parser_b.set_defaults(func=subcmd2_main)

        self.arg_parser = parser

        self.args = parser.parse_args()

    def run(self):
        self.run_arg_parser()

        if self.args.subcommand_name == 'subcmd1':
            subcmd1_main(self.args)
        elif self.args.subcommand_name == 'subcmd2':
            subcmd2_main(self.args)
        else:
            self.arg_parser.print_help()

def main():
    job = Job()
    job.run()

if __name__ == '__main__':
    main()