#!/usr/bin/env python

"""Writes MadGraph scripts to generate SM tt samples.

Choice of dynamic scale, PDF, and event selection are hard-coded.
"""

import argparse
import os


if __name__ == '__main__':
    
    argParser = argparse.ArgumentParser(epilog=__doc__)
    argParser.add_argument('directory', help='MadEvent directory')
    argParser.add_argument(
        '--mt', type=float, default=173.,
        help='Mass of top quark, GeV'
    )
    argParser.add_argument(
        '-n', type=int, default=5000000,
        help='Number of events to generate'
    )
    argParser.add_argument(
        '--split', type=int, default=100000,
        help='Number of events per LHE file'
    )
    argParser.add_argument(
        '-o', '--output', default='configs/ttbar',
        help='Prefix for the name of configuration files to be produced'
    )
    args = argParser.parse_args()
    
    
    outputDir = os.path.dirname(args.output)
    
    if outputDir and not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    
    for iClone in range(int(round(args.n / args.split))):
        
        jobName = os.path.basename('{}_{}'.format(args.output, iClone + 1))
        outFile = open('{}_{}.script'.format(args.output, iClone + 1), 'w')
        
        outFile.write('launch {} -n {}\n'.format(args.directory, jobName))
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
            outFile.write('{}{}\n'.format(indent, line))
        
        outFile.write('\n')
        outFile.close()
