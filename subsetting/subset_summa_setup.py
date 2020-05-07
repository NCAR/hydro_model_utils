#!/usr/bin/env python
"""Subset a SUMMA setup based on a list of hruIds

Parse an existing SUMMA file manager file and subset
the setup described by the file manager file based on 
a list of hruIds. Create a new SUMMA file manager file
that describes the subsetted setup.
"""

import argparse
from datetime import datetime
import os
import re
import sys
import xarray as xr

fm_keys = ['version', 'settings_path', 'forcing_path',
            'output_path', 'decisions', 'notused_1',
            'notused_2', 'notused_3', 'notused_4',
            'notused_5', 'output_control', 'notused_6',
            'notused_7', 'notused_8', 'local_attr',
            'local_param', 'basin_param', 'forcing_list',
            'initial_cond', 'param_trials', 'output_prefix']

otree = {'settings_path': 'settings',
         'forcing_path': 'forcings',
         'output_path': 'output'}

comment_sep = '!'

def parse_hruid_file(ifile):
    """Parse hruId file

    Parameters
    ----------
    ifile : str
        Pathname for text file with hruIds

    Returns
    -------
    list
        List of hruIds
    """
    with open(ifile) as f:
        ids = [int(x) for x in f]
    return ids

def fm_parse(ifile):
    """Parse SUMMA file manager file

    Parses a SUMMA file manager file to create a dictionary
    of its contents

    Parameters
    ----------
    ifile : str
        Pathname for SUMMA file manager file

    Returns
    -------
    dict
        file manager entries as key, value pairs
    dict
        file manager comments as key, value pairs
    str
        Content of original file
    """
    with open(ifile) as f:
        fm_txt = f.read()

    fm_values = []
    fm_comments = []
    for line in iter(fm_txt.splitlines()):
        m = re.match('^([^\\{}]*)\\{}(.*)$'.format(comment_sep, comment_sep), line)
        if m and m.group(1):  # The line contains a hash / comment
            fm_values.append(m.group(1).replace("'", ' ').strip())
            fm_comments.append(m.group(2))

    fm = dict(zip(fm_keys, fm_values))
    fm_comments = dict(zip(fm_keys, fm_comments))
 
    return fm, fm_comments, fm_txt
    
def process_command_line():
    """Parse the commandline

    Parameters
    ----------
    none

    Returns
    -------
    Namespace
        parsed commandline arguments
    """

    parser = argparse.ArgumentParser(description='Subset a SUMMA setup based on a list of hruIds.')
    parser.add_argument('filemanager',
                        help='filemanager file of the original SUMMA setup.')
    parser.add_argument('hruidfile', 
                        help='path of file with list of hruIds to subset.')
    parser.add_argument('opath',
                        help='directory where subsetted setup will be written.')

    args = parser.parse_args()

    return(args)

def fm_update(fm_org, opath, otree=otree):
    """Update file manager

    Parameters
    ----------
    fm_org : dict
        file manager to be updated

    opath : str
        output path for new SUMMA setup

    otree: dict (optional)
        directory tree for SUMMA setup

    Returns
    -------
    dict
        updated file manager
    """
    fm = fm_org.copy()
    for key in ('settings_path', 'forcing_path', 'output_path'):
        fm[key] = os.path.join(opath, otree[key])
        if fm[key][-1] != os.path.sep:
            fm[key] += os.path.sep

    return fm

def fm_write(fm, fm_comments, ofile, history=None):
    """Write the file manager to file

    Parameters
    ----------
    fm : dict
        Dictionary with filemanager entries
    fm_comments : dict
        Dictionary with filemanager comments
    ofile : str
        output file path
    
    Returns
    -------
    none
    """

    with open(ofile, 'w') as f:
        for key in fm:
            f.write("'{}' {}{}\n".format(fm[key], comment_sep, fm_comments[key]))
        if history:
            f.write("{} history: {}\n".format(comment_sep, history))

def otree_create(opath, otree=otree):
    """Create directory tree

    Parameters
    ----------
    opath : str
        Output path
    otree : dict (optional)
        Dictionary with output paths to create

    Returns
    -------
    none
    """
    opath = os.path.abspath(os.path.expanduser(opath))
    for sub in otree.values():
        os.makedirs(os.path.join(opath, sub), exist_ok=True)

def do_nothing(key, fm_org, fm, hru_ids):
    """Do nothing

    Parameters
    ----------
    key : str
        String to indicate entry in file manager
    fm_org : dict
        Dictionary with original file manager content
    fm : dict
        Dictionary with subsetted file manager content
    hru_ids : list
        List of HRU Ids


    Returns
    -------
    none
    """
    pass

def process_decisions(key, fm_org, fm, hru_ids):
    """Process the decisions file

    Copy the decisions file to the new location

    Parameters
    ----------
    key : str
        String to indicate entry in file manager
    fm_org : dict
        Dictionary with original file manager content
    fm : dict
        Dictionary with subsetted file manager content
    hru_ids : list
        List of HRU Ids


    Returns
    -------
    none
"""
    pass
    
processors = [do_nothing, do_nothing, do_nothing,
              do_nothing, process_decisions, do_nothing,
              do_nothing, do_nothing, do_nothing,
              do_nothing, do_nothing, do_nothing,
              do_nothing, do_nothing, do_nothing,
              do_nothing, do_nothing, do_nothing,
              do_nothing, do_nothing, do_nothing]

# processors = [do_nothing, do_nothing, do_nothing,
#               do_nothing, do_nothing, do_nothing,
#               do_nothing, process_output_control, do_nothing,
#               do_nothing, do_nothing, process_local_attr,
#               process_local_param, process_basin_param, process_forcing_list,
#               process_initial_cond, process_param_trials, process_output_prefix]

process = dict(zip(fm_keys, processors))

if __name__ == '__main__':

    # Parse input and prep for subsetting

    # process command line
    args = process_command_line()

    # make sure opath exists
    if not os.path.isdir(args.opath):
        sys.exit('{}: Output path does not exist'.format(args.opath))

    # parse the summa file manager
    fm_org, fm_org_comments, fm_org_txt = fm_parse(args.filemanager)
 
    # read the IDs to subset
    hru_ids = parse_hruid_file(args.hruidfile)

    # create a history string to be passed to all updated files
    history = '{}: {}\n'.format(datetime.now().strftime('%c'),
                                ' '.join(sys.argv))

    # go through the SUMMA filemanager and subset each of the files in turn.
    # write the subsetted files to opath using the directory structure 
    # defined by otree{} at the top of this script

    # update the filemanager for the new setup
    fm = fm_update(fm_org, args.opath)

    # write the file manager to file
    fm_write(fm, fm_org_comments, 
             os.path.join(args.opath, os.path.basename(args.filemanager)),
             history)

    # create the output directories according to otree
    otree_create(args.opath, otree)

    # process all the filemanager input
    for key in fm_keys:
        process[key](key, fm_org, fm, hru_ids)



