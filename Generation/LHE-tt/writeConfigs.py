#!/usr/bin/env python

"""Writes MadGraph scripts to generate SM tt samples.

Choice of dynamic scale, PDF, and event selection are hard-coded.
"""

import argparse
import os


if __name__ == '__main__':
    
    arg_parser = argparse.ArgumentParser(epilog=__doc__)
    arg_parser.add_argument('directory', help='MadEvent directory')
    arg_parser.add_argument(
        '--mt', type=float, default=173.,
        help='Mass of top quark, GeV'
    )
    arg_parser.add_argument(
        '-n', type=int, default=5000000,
        help='Number of events to generate'
    )
    arg_parser.add_argument(
        '--split', type=int, default=100000,
        help='Number of events per LHE file'
    )
    arg_parser.add_argument(
        '-o', '--output', default='configs/ttbar',
        help='Prefix for the name of configuration files to be produced'
    )
    args = arg_parser.parse_args()
    
    
    output_dir = os.path.dirname(args.output)
    
    if output_dir:
        try:
            os.makedirs(output_dir)
        except FileExistsError:
            pass
    
    
    for iclone in range(int(round(args.n / args.split))):
        
        job_name = os.path.basename('{}_{}'.format(args.output, iclone + 1))
        outfile = open('{}_{}.script'.format(args.output, iclone + 1), 'w')
        
        outfile.write('launch {} -n {}\n'.format(args.directory, job_name))
        indent = ' ' * 2
        
        for line in [
            '',
            'done',
            '',
            'set run_card nevents {}'.format(args.split),
            'set run_card pdlabel lhapdf',
            'set run_card lhaid 90400',
            'set run_card sys_pdf PDF4LHC15_nlo_30_pdfas',
            'set run_card dynamical_scale_choice 0',
            '',
            'set run_card xptl 30',
            '',
            'set param_card mt  {}'.format(args.mt),
            'set param_card ymt {}'.format(args.mt),
            '',
            '',
            'done',
        ]:
            outfile.write('{}{}\n'.format(indent, line))
        
        outfile.write('\n')
        outfile.close()
