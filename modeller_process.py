import pypdb as pb
import urllib.request
import urllib.error
import xmltodict
from json import loads, dumps
import os
import time
from Bio import PDB
import pypdb as pb

from modeller import *
from modeller.automodel import *


class SelectChain(PDB.Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0

    def accept_chain(self, chain):
        return chain.get_id() == self.chain



def get_chain_pdb(id, chain):
    cache_path = './cache/' + id + chain + '.pdb'
    if os.path.isfile(cache_path):
        return cache_path
    parser = PDB.PDBParser()
    writer = PDB.PDBIO()
    writer.set_structure(parser.get_structure(id, get_pdb(id)))
    writer.save(cache_path, select=SelectChain(chain))
    return cache_path


def get_pdb(id):
    if len(id) == 5:
        return get_chain_pdb(id[:4], id[4])
    cache_path = './cache/' + id + '.pdb'
    if os.path.isfile(cache_path):
        return cache_path
    if len(id) != 4:
        print("BAD ID:", id)
    with open(cache_path, 'a') as file:
        file.write(pb.get_pdb_file(id))
    return cache_path

from subprocess import Popen, PIPE

def tm_align(id1, id2):
    cache_path = './cache/' + id1 + '_' + id2 + '.ali'
    if os.path.isfile(cache_path):
        with open(cache_path, 'r') as file:
            return file.read()
    process = Popen(["./tmalign_folder/TMalign", get_pdb(id1), get_pdb(id2)],
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    result = stdout.decode("utf-8").split('\n')
    if len(result) == 27:
        seq1 = result[22]
        seq2 = result[24]
        rmsd = result[16].split("RMSD=")[1].split(',')[0].split(" ")[-1]
        result = []
        result.append("C;" + rmsd)
        result.append(">P1;" + id1)
        result.append("sequence:" + id1 + ":.:.:.:.::::")
        result.append(seq1 + '*')
        result.append(">P1;" + id2)
        result.append("structure:" + id2 + ":.:.:.:.::::")
        result.append(seq2 + '*')
        result = '\n'.join(result)
        with open(cache_path, 'a') as file:
            file.write(result)
    else:
        result = ""
    return result



def get_rmsd(id1, id2):
    result = tm_align(id1, id2)
    return result.split('\n')[0][2:-1]



def generate_model(target, template):
    if target == template:
        return
    log.verbose()
    env = environ()
    env.io.atom_files_directory = ['.', './cache']
    allow_alternates=True
    a = automodel(env,
                  alnfile  = './cache/' + target + '_' + template + '.ali',
                  knowns   = template,
                  sequence = target)
    a.starting_model= 1
    a.ending_model  = 1
    a.make()

    os.rename(target + '.B99990001.pdb', './cache/' + target + '.model' + template + '.pdb')
    rmsd = get_rmsd(target, target + '.model' + template)

    directory = './cache/models/' + target + '/' + rmsd + '/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    os.rename('./cache/' + target + '.model' + template + '.pdb', directory + template + '.pdb')
    os.remove(target + '.D00000001')
    os.remove(target + '.ini')
    os.remove(target + '.rsr')
    os.remove(target + '.sch')
    os.remove(target + '.V99990001')

import sys
if __name__ == "__main__":
    generate_model(sys.argv[1], sys.argv[2])
