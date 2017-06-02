#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################################
#                                                                                                           #
#           re_process_hrc_data.py: a control script to run reprocess csh scripts: HRC S Version            #
#                                                                                                           #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                                       #
#                                                                                                           #
#           Last Update: Mar 09, 2017                                                                       #
#                                                                                                           #
#############################################################################################################

import sys
import os
import string
import re
import math
import unittest
import time
import numpy
import astropy.io.fits  as pyfits
from datetime import datetime
#
#--- from ska
#
from Ska.Shell import getenv, bash

ascdsenv = getenv('source /home/ascds/.ascrc -r release; source /home/mta/bin/reset_param ', shell='tcsh')
ciaoenv  = getenv('source /soft/ciao/bin/ciao.sh')

#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project9/Scripts/Hrc_S/house_keeping/dir_list'

f    = open(path, 'r')
data = [line.strip() for line in f.readlines()]
f.close()

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec "%s = %s" %(var, line)
#
#--- append path to a private folders
#
sys.path.append(mta_dir)
sys.path.append(bin_dir)

import mta_common_functions     as mcf
import convertTimeFormat        as tcnv
import OcatSQL                  as sql
from   OcatSQL                  import OcatDB

#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

#-----------------------------------------------------------------------------------------
#-- run_process: a control script to run reprocess csh scripts                         ---
#-----------------------------------------------------------------------------------------

def run_process():
    """
    a control script to run reprocess csh scripts 
    input:  none
    output: hrc_i_list  --- a list of hrc i obsids which need to be re-processed
            hrc_s_list  --- a list of hrc s obsids which need to be re-processed
            <data_dir>/<obsid>    --- re-processed data direcotry
    """

    [hrc_i, hrc_s] = find_un_processed_data()

#    print "HRC I : " + str(hrc_i)
    print "HRC S : " + str(hrc_s)

#    fo  = open('hrc_i_list', 'w')
#    for ent in hrc_i:
#        fo.write(str(ent))
#        fo.write('\n')
#    fo.close()
#
    fo  = open('hrc_s_list', 'w')
    for ent in hrc_s:
        fo.write(str(ent))
        fo.write('\n')
    fo.close()

#    cmd = 'csh -f ' + bin_dir + 'repro_all.csh hrc_i_list'
#    run_ciao(cmd)

    cmd = 'csh -f ' + bin_dir + 'repro_all_S.csh hrc_s_list'
    run_ciao(cmd)

    for obsid in hrc_s:
        cmd = 'chgrp -R hat ' +  data_dir + '/' + obsid
        os.system(cmd)


#-----------------------------------------------------------------------------------------
#-- find_un_processed_data: find hrc obsids which need to be reprocessed                --
#-----------------------------------------------------------------------------------------

def find_un_processed_data():
    """
    find hrc obsids which need to be reprocessed
    input: none
    output: uhrc_i/uhrc_s   --- lists of obsids of hrc i and hrc s
    """
#
#--- extract all hrc obsid listed in database
#
    infile = '/data/mta4/obs_ss/sot_ocat.out'
    data   = read_data(infile)
    hrc_i  = []
    hrc_s  = []
    for ent in data:
        atemp = re.split('\^', ent)
        mc1 = re.search('HRC-I', atemp[12])
        mc2 = re.search('HRC-S', atemp[12])
        if mc1 is not None:
            atemp = re.split('\^\s+', ent)
            atemp[0].strip()
            try:
                val   = int(float(atemp[1]))
            except:
                continue
            hrc_i.append(val)

        elif mc2 is not None:
            atemp = re.split('\^\s+', ent)
            atemp[0].strip()
            try:
                val   = int(float(atemp[1]))
            except:
                continue
            hrc_s.append(val)

        else:
            continue

    uhrc_i = clean_the_list(hrc_i, 'i')
    uhrc_s = clean_the_list(hrc_s, 's')
    
    return [uhrc_i, uhrc_s]

#-----------------------------------------------------------------------------------------
#-- find_processed_data: find the hrc obsids which are already re-processed             --
#-----------------------------------------------------------------------------------------

def find_processed_data(inst):
    """
    find the hrc obsids which are already re-processed
    input:  inst    --- instrument designation: "i" or "s"
    output: out     --- a list of obsids
    """
    cmd  = 'ls -d ' + data_dir + '/* > ' + zspace
    os.system(cmd)
    data = read_data(zspace, remove=1)

    out = []
    for ent in data:
        atemp = re.split('\/', ent)
        try:
            val   = int(float(atemp[-1]))
        except:
            continue
        if mcf.chkNumeric(val):
            out.append(val)
#
#--- remove duplicate
#
    oset = set(out)
    out  = list(oset)

    return out

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

def clean_the_list(current, inst):
    """
    select out obsids which need re-process
    input:  current --- a list of all hrc obsid found in database
            inst    --- instrument designation; "i" or "s"
    output: good    --- a list of obsids to be re-processed
            <house_keeping>/cancelled_list  --- a list of observations canceleld or discarded
    """
#
#--- read the past cancelled obsids
#
    rfile  = house_keeping + 'cancelled_list'
    remove = set(read_data(rfile))
#
#--- find obsids already re-processed
#
    phrc = find_processed_data(inst)
    uhrc = set(current) - set(phrc)
    uhrc = list(uhrc - remove)
#
#--- select out obsids which need to be reprocessed
#
    good = []
    bad  = []
    for obsid in uhrc:
        status = find_status(obsid)
        if status == 'archived':
            good.append(obsid)
        elif status == 'observed':
            pass
        else:
            bad.append(obsid)
#
#--- update cancelled_list if there are new cancelled observations
#
    if len(bad) > 0:
        fo = open(rfile, 'a')
        for ent in  bad:
            fo.write(str(ent))
            fo.write('\n')

        fo.close()

    return good

#-----------------------------------------------------------------------------------------
#-- find_status: find the status of the observations from the database                  --
#-----------------------------------------------------------------------------------------

def find_status(obsid):
    """
    find the status of the observations from the database
    input:  obsid   --- obsid
    output: status  --- status
    """
    try:
        dbase  = OcatDB(obsid)
        status = dbase.origValue('status')
        return status
    except:
        return ""

#-----------------------------------------------------------------------------------------
#-- read_data: read data file                                                           --
#-----------------------------------------------------------------------------------------

def read_data(infile, remove=0):

    f    = open(infile, 'r')
    data = [line.strip() for line in f.readlines()]
    f.close()

    if remove == 1:
        mcf.rm_file(infile)

    return data

#-----------------------------------------------------------------------------------------
#-- run_ciao: running ciao comannds                                                    ---
#-----------------------------------------------------------------------------------------

def run_ciao(cmd, clean =0):
    """
    run the command in ciao environment
    input:  cmd --- command line
    clean   --- if 1, it also resets parameters default: 0
    output: command results
    """
    if clean == 1:
        acmd = '/usr/bin/env PERL5LIB=""  source /home/mta/bin/reset_param ;' + cmd
    else:
        acmd = '/usr/bin/env PERL5LIB="" LD_LIBRARY_PATH=""   ' + cmd
    
    try:
        bash(acmd, env=ciaoenv)
    except:
        try:
            bash(acmd, env=ciaoenv)
        except:
            pass

#-----------------------------------------------------------------------------------------

if __name__ == '__main__':

    run_process()
