#!/usr/bin/env python

"""Writes MadGraph scripts to generate signal samples.

Choice of dynamic scale, PDF, event selection, and masses of the heavy
Higgs bosons are hard-coded.
"""

import argparse
import os


if __name__ == '__main__':
    
    arg_parser = argparse.ArgumentParser(epilog=__doc__)
    arg_parser.add_argument('state', help='CP state, A or H')
    arg_parser.add_argument('directory', help='MadEvent directory')
    arg_parser.add_argument(
        '-w', '--width', type=float, default=0.1,
        help='Relative total width'
    )
    arg_parser.add_argument(
        '-n', type=int, default=1000000,
        help='Number of events to generate'
    )
    arg_parser.add_argument(
        '--split', type=int, default=100000,
        help='Number of events per LHE file'
    )
    arg_parser.add_argument(
        '-o', '--output', default='configs',
        help='Name for directory to store configurations'
    )
    args = arg_parser.parse_args()
    
    if args.state not in ['A', 'H']:
        raise RuntimeError('Unsupported CP state "{}".'.format(args.state))
    
    if args.output:
        try:
            os.makedirs(args.output)
        except FileExistsError:
            pass
    
    
    if args.state == 'A':
        param_mass = 'MA0'
        param_width = 'A0Width'
        param_pdgid = 6000046
    else:
        param_mass = 'MH0'
        param_width = 'H0Width'
        param_pdgid = 6000045
    
    
    for mass in [400, 500, 600, 700, 1000]:
        width = mass * ars.width
            
        for iclone in range(int(round(args.n / args.split))):
            
            job_name = 'm{:g}_w{:g}_{}'.format(mass, width, iclone + 1)
            outfile = open(os.path.join(args.output, job_name + '.script'), 'w')
            
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
                'set param_card {}  {}'.format(param_mass, mass),
                'set param_card {}  {}'.format(param_width, width),
                'set param_card DECAY {}  {}'.format(param_pdgid, width),
                '',
                'done',
            ]:
                outfile.write('{}{}\n'.format(indent, line))
            
            outfile.write('\n')
            outfile.close()
