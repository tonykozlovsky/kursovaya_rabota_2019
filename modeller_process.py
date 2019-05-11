import pypdb as pb
import urllib.request
import urllib.error
import xmltodict
from json import loads, dumps
import os
import time
from Bio import PDB
import pypdb as pb
import sys

from modeller import *
from modeller.automodel import *

def console(args):
    return Popen(args.split(), stdout=PIPE).communicate()[0].decode("utf-8")

cwd = "kek"

def get_chain_pdb(id, chain):
    global cwd
    cache_path = cwd + '/cache/data/' + id + chain + '.pdb'
    if os.path.isfile(cache_path):
        return cache_path
    sys.stderr.write("CACHE MISS WHEN REQUIRED: " + id + chain)
    sys.exit(228)


def get_pdb(id):
    global cwd
    if len(id) > 4:
        return get_chain_pdb(id[:4], id[4:])
    cache_path = cwd + '/cache/data/' + id + '.pdb'
    if os.path.isfile(cache_path):
        return cache_path
    sys.stderr.write("CACHE MISS WHEN REQUIRED: " + id)
    sys.exit(228)


from subprocess import Popen, PIPE


def get_rmsd(id1, id2):
    global cwd
    process = Popen([cwd + "/tmalign_folder/TMalign", get_pdb(id1), get_pdb(id2)],
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    result = stdout.decode("utf-8").split('\n')
    if len(result) == 27:
        return result[16].split("RMSD=")[1].split(',')[0].split(" ")[-1]
    sys.stderr.write("CANNOT GET RMSD: " + id1 + " " + id2)
    sys.exit(228)


def generate_model(target, template, arg_cwd):
    global cwd
    cwd = arg_cwd
    if target == template:
        return
    alnfile = cwd + '/cache/data/' + target + '_' + template + '.ali'
    if not os.path.isfile(alnfile):
        sys.stderr.write("NO ALN FILE: " + alnfile)
        sys.exit(228)
    log.verbose()
    env = environ()
    env.io.atom_files_directory = ['.', cwd + '/cache/data']
    a = automodel(env,
                  alnfile  = alnfile,
                  knowns   = template,
                  sequence = target)
    a.starting_model= 1
    a.ending_model  = 1
    a.make()

    os.rename(target + '.B99990001.pdb', cwd + '/cache/data/' + target + '.model' + template + '.pdb')
    rmsd = get_rmsd(target, target + '.model' + template)

    directory = cwd + '/cache/models/' + target + '/' + rmsd + '/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    os.rename(cwd + '/cache/data/' + target + '.model' + template + '.pdb', directory + template + '.pdb')
    console("pigz --best -f -q " +  directory + template + '.pdb')
    os.remove(target + '.D00000001')
    os.remove(target + '.ini')
    os.remove(target + '.rsr')
    os.remove(target + '.sch')
    os.remove(target + '.V99990001')


if __name__ == "__main__":
    generate_model(sys.argv[1], sys.argv[2], sys.argv[3])
