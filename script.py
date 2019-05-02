import pypdb as pb
import urllib.request
import urllib.error
import xmltodict
from json import loads, dumps
import os
import time
from Bio import PDB


def remove_at_sign(kk):
    tagged_keys = [thing for thing in kk.keys() if thing.startswith('@')]
    for tag_key in tagged_keys:
        kk[tag_key[1:]] = kk.pop(tag_key)
    return kk


def to_dict(odict):
    out = loads(dumps(odict))
    return out


def get_info(pdb_id, url_root='http://www.rcsb.org/pdb/rest/describeMol?structureId='):
    url = url_root + pdb_id
    req = urllib.request.Request(url)
    f = urllib.request.urlopen(req)
    result = f.read()
    assert result
    out = xmltodict.parse(result, process_namespaces=True)
    return out

def query(q_type, arn_name, value):
    url = "http://www.rcsb.org/pdb/rest/" + q_type + "?" + arn_name + "=" + str(value)
    req = urllib.request.Request(url)
    f = urllib.request.urlopen(req)
    result = f.read()
    out = xmltodict.parse(result, process_namespaces=True)
    return to_dict(out)


def get_clusterr_domains(pdb_id):
    out = get_info(pdb_id, url_root =
        'http://www.rcsb.org/pdb/rest/representatives?structureId=')
    out = to_dict(out)
    return remove_at_sign(out['representatives'])

import xml

def get_candidates(id, chains):
    cache_path = './cache/' + id + chains[0] + '.candidates'

    if os.path.isfile(cache_path):
        with open(cache_path, 'r') as file:
            return file.read()

    try:
        response = pb.get_seq_cluster(id + '.' + chains[0])['pdbChain']
    except xml.parsers.expat.ExpatError:
        response = []
    result = []

    if type(response) != type([]):
        response = [response]

    for el in response:
        name = el['@name'].split('.')
        result.append(name[0] + ' ' + name[1] + ' ' + el['@rank'])

    result = '|'.join(result)

    with open(cache_path, 'a') as file:
        file.write(result)

    return result


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
    if not len(result) < 25:
        seq1 = result[18]
        seq2 = result[20]
        rmsd = result[12].split("RMSD=")[1].split(',')[0].split(" ")[-1]
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
    #else:
        #print(result)
    return result


class SelectChain(PDB.Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_residue(self, residue):
        return len(str(residue).split("het= ")) > 1

    def accept_chain(self, chain):
        return chain.get_id() == self.chain

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


def get_chain_pdb(id, chain):
    cache_path = './cache/' + id + chain + '.pdb'
    if os.path.isfile(cache_path):
        return cache_path
    parser = PDB.PDBParser()
    writer = PDB.PDBIO()
    writer.set_structure(parser.get_structure(id, get_pdb(id)))
    writer.save(cache_path, select=SelectChain(chain))
    return cache_path


def get_chains(id):
    cache_path = './cache/' + id + '.chains'

    if os.path.isfile(cache_path):
        with open(cache_path, 'r') as file:
            return file.read()

    chains = []
    response = pb.get_entity_info(id)['Entity']['Chain']
    if type(response) != type([]):
        response = [response]

    for el in response:
        chains.append(el['@id'])

    chains = "".join(chains)

    with open(cache_path, 'a') as file:
        file.write(chains)

    return chains


def get_rmsd(id1, id2):
    result = tm_align(id1, id2)
    return result.split('\n')[0][2:-1]


import threading

global_start_time = 0
sum_time = 0
models = 0
bad_models = 0
cache_models = 0

def write_error(target, template):
    with open('./error_candidates.txt', 'a') as file:
        file.write(target + " " + template + "\n")

def process_id(id, all_models):
    global global_start_time
    global sum_time
    global threads
    global models
    global bad_models
    global cache_models
    get_pdb(id)
    candidates = get_candidates(id, get_chains(id)).split("|")
    target = id + get_chains(id)[0]
    for candidate in candidates:
        start_time = time.time()
        candidate = candidate.split(" ")
        if len(candidate) < 2:
            continue
        template = candidate[0] + candidate[1]
        try:
            if target == template:
                continue
            if (target, template) in all_models:
                cache_models += 1
                continue
            get_chain_pdb(candidate[0], candidate[1])
            alignment = tm_align(target, template)
            if len(alignment) < 25:
                write_error(target, template)
                print("ALIGNMENT ERROR")
                bad_models += 1
            else:
                process = Popen(["python3", "./modeller_process.py",
                                 target, template], stdout=PIPE, stderr=PIPE)
                stdout, stderr = process.communicate()
                if len(stderr) > 0:
                    write_error(target, template)
                    print("MODELLER ERROR")
                    bad_models += 1
                else:
                    models += 1
                    sum_time += time.time() - start_time
                    print("new model:", target, template, "   ",
                          "cur_time:", time.time() - start_time,
                          "   ", "models:", str(models) + "/" + str(bad_models) + "/"
                          + str(cache_models),
                          "   ", "time per model:", (time.time() - global_start_time) / models)
        except Exception as e:
            # write_error(target, template)
            # print("EXCEPTION", e)
            bad_models += 1
            continue

def clean():
    stdout = os.popen("find ./ -type f").read().split("\n")
    stdout.pop()
    for file in stdout:
        if file.endswith(".rsr") or file.endswith(".ini") or file.endswith(".D00000001") or \
                file.endswith(".sch"):
            os.popen("rm " + file)
            print("removed: " + file)
        if file.startswith("./cache/") and len(file.split("/")) == 3:
            os.popen("rm " + file)
            print("removed: " + file)
    return

def get_error_models():
    res = []
    with open("./error_candidates.txt", 'r') as file:
        res = file.read().split("\n")
        res.pop()
    error_candidates = []
    for line in res:
        ec = line.split(" ")
        error_candidates.append((ec[0], ec[1]))
    return error_candidates

def get_all_generated_models():
    result = set()
    stdout = os.popen("find ./cache/models -type f").read().split("\n")
    stdout.pop()
    for model in stdout:
        model = model.split("/")
        result.add((model[3], model[5].split(".")[0]))
    for pair in get_error_models():
        result.add(pair)
    return result

def main():
    global global_start_time
    #clean()
    all_models = get_all_generated_models()
    print("All models:", len(all_models))
    time.sleep(5)
    threads_pool = set()
    global_start_time = time.time()
    with open('./data/sids.txt', 'r') as file:
        for id in file.read().split('\n'):
            while True:
                new_threads_pool = set()
                for t in threads_pool:
                    if t.is_alive():
                        new_threads_pool.add(t)
                threads_pool = new_threads_pool
                if len(new_threads_pool) < 64:
                    break
                time.sleep(0.01)
            # print("NEW THREAD:", len(threads_pool) + 1)
            t = threading.Thread(target=process_id, args=(id,all_models,))
            t.start()
            threads_pool.add(t)

    return

    #t = time.time()
    #c = time.clock()
    #do_alignment('5GW9.A', '5ZO6.X')
    #print("Time: " + str(time.time() - t))
    #print("CPU: " + str(time.clock() - c))

    #generate_for_chain('5ZOH')
    #print(get_clusterr_domains('5ZOH'))


    #print(rmsd.kabsch(
    #    rmsd.get_coordinates_pdb(id1 + '.pdb')[1],
    #    rmsd.get_coordinates_pdb(id2 + '.pdb')[1]))


    #print(query('sequenceCluster', 'structureId', '5ZOH.A'))
    #print(query('representatives', 'structureId', '5ZOH.A'))
    #print(query('representativeDomains', 'structureId', '5ZOH.A'))

if __name__ == "__main__":
    main()
