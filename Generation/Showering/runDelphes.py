#!/usr/bin/env python

"""Runs Delphes simulation starting from LHE files.

Arguments to this script is a list of directories with results of
indivudal MadEvent runs (i.e. each of these directories contains file
"unweighted_events.lhe.gz").  For each of the thus provided LHE file the
script runs Pythia and Delphes using specified configuration files.
The processing is done in parallel.  Produced Delphes files are placed
in given output directory (defaults to the current one).  They are named
after the base names of the input directories, and user can also specify
an optional prefix.  For each job a log file is written and saved in
subdirectory "logs" in the output directory.
"""

import argparse
from datetime import datetime
import os
import queue
import re
import subprocess
import sys
import threading
from uuid import uuid4


# Global lock for printing
printLock = threading.Lock()


def worker(inputs, delphesConfigFile, pythiaConfigTemplate, prefix, outDir, logDir, tmpDir):
    
    while True:
        try:
            lheFile = inputs.get_nowait()
        except queue.Empty:
            break
        
        
        # Construct base name for the output file and log file
        name = prefix + os.path.basename(os.path.dirname(lheFile))
        outFilePath = os.path.join(outDir, name + '.root')
        
        
        # Create log file
        logFile = open(os.path.join(logDir, name + '.log'), 'w')
        logFile.write('Job to produce file "{}" started at {}.\n\n'.format(
            outFilePath, datetime.now()
        ))
        
        
        # Make sure the input file actually exists
        if not os.path.exists(lheFile) or not os.path.isfile(lheFile):
            logFile.write('File "{}" does not exist. Job aborted.\n'.format(lheFile))
            logFile.close()
            
            with printLock:
                print('File "{}" does not exist.'.format(lheFile), file=sys.stderr)
            
            inputs.task_done()
            return
            
        
        
        tmpFileNames = []
        
        # Uncompress the LHE file into a temporary copy
        unzippedLHEName = os.path.join(tmpDir, 'events_{}.lhe'.format(name))
        unzippedLHE = open(unzippedLHEName, 'w')
        subprocess.check_call(
            ['gzip', '-cd', lheFile], stdout=unzippedLHE
        )
        unzippedLHE.close()
        
        logFile.write(
            'Input LHE file "{}" unzipped to file "{}".\n\n'.format(
                lheFile, unzippedLHEName
            )
        )
        tmpFileNames.append(unzippedLHEName)
        
        
        # Create a temporary file with Pythia configuration
        pythiaConfigFileName = os.path.join(tmpDir, 'pythiaConfig_{}.cmnd'.format(name))
        pythiaConfigFile = open(pythiaConfigFileName, 'w')
        pythiaConfigFile.write(pythiaConfigTemplate)
        pythiaConfigFile.write('Beams:LHEF = {}\n'.format(unzippedLHEName))
        pythiaConfigFile.close()
        
        logFile.write('Temporary configuration file "{}" created.\n\n'.format(
            pythiaConfigFileName
        ))
        tmpFileNames.append(pythiaConfigFileName)
        
        
        logFile.flush()
        
        # Run Delphes
        logFile.write('Starting Delphes.\n')
        
        command = ['DelphesPythia8', delphesConfigFile, pythiaConfigFileName, outFilePath]
        
        try:
            subprocess.check_call(command, stdout=logFile, stderr=subprocess.STDOUT)
        
        except subprocess.CalledProcessError as error:
            logFile.write('Command "{}" terminated with error code {}. Job aborted.\n'.format(
                ' '.join(command), error.returncode
            ))
            logFile.close()
            
            with printLock:
                print(
                    'Delphes terminated with an error when processing file "{}".'.format(lheFile),
                    file=sys.stderr
                )
            
            inputs.task_done()
            return
        
        logFile.write('\nDelphes run is complete.\n\n')
        
        
        # Clean up temporary files
        for f in tmpFileNames:
            os.remove(f)
            logFile.write('Temporary file "{}" deleted.\n'.format(f))
        
        
        with printLock:
            print('Finished processing file "{}".'.format(lheFile))
        
        logFile.write('\nEverything done.\n')
        logFile.close()
        inputs.task_done()


if __name__ == '__main__':
    
    argParser = argparse.ArgumentParser(epilog=__doc__)
    argParser.add_argument(
        'inputDirs', nargs='+',
        help='Directories that store results of individual MadEvent runs'
    )
    argParser.add_argument(
        '-o', '--output', default='.',
        help='Directory to store produced files'
    )
    argParser.add_argument(
        '-p', '--prefix', default='',
        help='Prefix to be added to names of output files'
    )
    argParser.add_argument(
        '-n', '--num-parallel', type=int, dest='numParallel', default=16,
        help='Number of jobs to run in parallel'
    )
    argParser.add_argument(
        '--pythia-config', dest='pythiaConfig', default='pythiaConfig.cmnd',
        help='Configuration file for Pythia'
    )
    argParser.add_argument(
        '--delphes-config', dest='delphesConfig', default='delphes_card.tcl',
        help='TCL configuration file for Delphes'
    )
    args = argParser.parse_args()
    
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    logDir = os.path.join(args.output, 'logs')
    
    if not os.path.exists(logDir):
        os.makedirs(logDir)
    
    for path in [args.pythiaConfig, args.delphesConfig]:
        if not os.path.exists(path) or not os.path.isfile(path):
            raise RuntimeError('File "{}" does not exist.'.format(path))
    
    
    startTime = datetime.now()
    
    inputs = queue.Queue()
    
    for d in args.inputDirs:
        inputs.put(os.path.join(d, 'unweighted_events.lhe.gz'))
    
    
    # Create a Pythia configuration but without the name of the input
    # LHE file, which will be added in each job
    with open(args.pythiaConfig) as f:
        pythiaConfigTemplate = f.readlines()
    
    for i in range(len(pythiaConfigTemplate)):
        if 'Beams:LHEF' in pythiaConfigTemplate[i]:
            del pythiaConfigTemplate[i]
            break
    
    pythiaConfigTemplate = ''.join(pythiaConfigTemplate)
    
    
    # Create a directory to store temporary files
    tmpDir = uuid4().hex
    os.makedirs(tmpDir)
    
    
    # Run generation using a thread pool
    threads = []
    
    for i in range(min(args.numParallel, inputs.qsize())):
        t = threading.Thread(target=worker, args=(
            inputs, args.delphesConfig, pythiaConfigTemplate,
            args.prefix, args.output, logDir, tmpDir
        ))
        threads.append(t)
        t.start()
    
    for t in threads:
        t.join()
    
    
    elapsedTime = datetime.now() - startTime
    
    try:
        os.rmdir(tmpDir)
        print('Done. Total elapsed time: {}.'.format(elapsedTime))
    except OSError:
        print(
            'Processing done but temporary directory "{}" is not emptry. '
            'Some jobs must have failed.'.format(tmpDir)
        )
        print('Total elapsed time: {}.'.format(elapsedTime))
