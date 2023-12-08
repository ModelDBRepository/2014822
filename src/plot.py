from src.plots.general import plot_by_datestring
import argparse
import os
from glob import glob
from pathlib import Path

class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        return argparse.HelpFormatter._split_lines(self, text, width)

def parse_args():
    parser = argparse.ArgumentParser(description='Plot runs',
                                    formatter_class=SmartFormatter)

    parser.add_argument("-p", "--protocol",
                            choices=[
                                "all",
                                ],
                            default="all",
                            help="What to plot"
                            )
    parser.add_argument("-i", "--id", help="Which experiment to run from latest", default="0")
        
    args = parser.parse_args()

    Path(f"figures").mkdir(parents=True, exist_ok=True)
    print(f"Plotting {args.protocol} experiment {args.id} ...")
    print("------------------------")

    return args

def all(args):
    protocol = "all"
    experiment_id = int(args.id)
    allparams = ['ltp', 'ltd']
    args.suffix = ""
    args.protocol = "all"

    print(f"Reading debug {protocol} experiment id: {experiment_id}...")
    files = list(filter(os.path.isdir, glob(f"debug/{protocol}/*")))

    files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    files = [os.path.basename(x) for x in files]

    datetime = list(dict.fromkeys([".".join(x.split('.')[:-1]) for x in files]))[experiment_id]
    filenames = [f for f in files if datetime in f]

    plot_by_datestring(filenames, args.suffix, args.protocol, allparams)


def main():
    if "SHOW_PLOT" not in os.environ:
        os.environ['SHOW_PLOT'] = 'FALSE'
    if "MILEDIDEBUG" not in os.environ:
        os.environ['MILEDIDEBUG'] = 'FALSE'
    

    args = parse_args()
    globals()[args.protocol](args)


if __name__ == "__main__":
    main()