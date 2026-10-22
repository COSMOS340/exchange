"""Explicit commands for the independent full-epithelial task."""
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'results/code'))


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('verify-input')
    sub.add_parser('prepare')
    fit = sub.add_parser('fit')
    fit.add_argument('--method', choices=['Monocle2', 'Slingshot'], required=True)
    pack = sub.add_parser('pack-return')
    pack.add_argument('--label',required=True)
    args = parser.parse_args()
    if args.command == 'verify-input':
        from pipeline_common import verify_data
        verify_data()
    elif args.command == 'prepare':
        from prepare_all import main as prepare
        prepare()
    elif args.command == 'fit':
        from run_fit_all import main as fit_main
        return fit_main(args.method)
    elif args.command == 'pack-return':
        from pack_return_all import pack_return
        pack_return(args.label)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
