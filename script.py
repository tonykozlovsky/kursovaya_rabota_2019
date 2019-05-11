import pypdb as pb
import urllib.request
import urllib.error
import xmltodict
from json import loads, dumps
import os
import time
from Bio import PDB
import pypdb as pb
from subprocess import Popen, PIPE
from threading import Thread


from modeller import *
from modeller.automodel import *

cwd = os.getcwd()


def console(args):
    return Popen(args.split(), stdout=PIPE).communicate()[0].decode("utf-8")


class SelectChain(PDB.Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0

    def accept_chain(self, chain):
        return chain.get_id() == self.chain

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


def get_chain_pdb(id, chain, cache):
    cache_path = './cache/data/' + id + chain + '.pdb'
    if os.path.isfile(cache_path):
        return cache_path
    if cache:
        print("CACHE MISS WHEN REQUIRED: ", id, chain)
    parser = PDB.PDBParser()
    writer = PDB.PDBIO()
    writer.set_structure(parser.get_structure(id, get_pdb(id, cache)))
    writer.save(cache_path, select=SelectChain(chain))
    return cache_path


def get_chains(id):
    cache_path = './cache/data/' + id + '.chains'

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

    with open(cache_path, 'w') as file:
        file.write(chains)

    return chains

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
    cache_path = './cache/data/' + id + chains[0] + '.candidates'

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

    with open(cache_path, 'w') as file:
        file.write(result)

    return result


def get_pdb(id, cache):
    if len(id) > 4:
        return get_chain_pdb(id[:4], id[4:], cache)
    cache_path = './cache/data/' + id + '.pdb'
    if os.path.isfile(cache_path):
        return cache_path
    if len(id) != 4:
        print("BAD ID:", id)
    if cache:
        print("CACHE MISS WHEN REQUIRED: ", id)
    with open(cache_path + ".gz", 'wb') as file:
        file.write(pb.get_pdb_file(id, compression=True))
    console("unpigz -f -q " + cache_path + ".gz")
    return cache_path


def tm_align(id1, id2, cache):
    cache_path = './cache/data/' + id1 + '_' + id2 + '.ali'
    if os.path.isfile(cache_path):
        with open(cache_path, 'r') as file:
            return file.read()

    process = Popen(["./tmalign_folder/TMalign", get_pdb(id1, cache), get_pdb(id2, cache)],
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
        with open(cache_path, 'w') as file:
            file.write(result)
    else:
        result = ""
    return result

def clean():
    stdout = console("find ./ -type f").split("\n")
    stdout.pop()
    for file in stdout:
        if file.endswith(".rsr") or file.endswith(".ini") or file.endswith(".D00000001") or \
                file.endswith(".sch") or file.endswith(".V99990001"):
            console("rm " + file)
            print("removed: " + file)
        #if file.endswith(".gz.gz"):
        #    console("rm " + file)
        #    print("removed: " + file)
        #if file.endswith(".pdb") and file.startswith("./cache/models/"):
        #    console("rm " + file)
        #    print("removed: " + file)
    stdout = console("find ./cache/data -type f").split("\n")
    stdout.pop()
    for file in stdout:
        console("rm " + file)
        print("removed: " + file)
        pass
    if os.path.exists("./scripts"):
        console("rm -r ./scripts")
    return

def get_all_generated_models():
    result = set()
    stdout = console("find ./cache/models -type f").split("\n")
    stdout.pop()
    for model in stdout:
        model = model.split("/")
        result.add((model[3], model[5].split(".")[0]))
    return result

from mpipe import Pipeline, UnorderedStage

def run_modeller(arg):
    start = time.time()
    target, template = arg.split()
    alignment = tm_align(target, template, True)
    if len(alignment) == 0:
        print("ALIGNMENT ERROR")
    else:
        directory = './scripts/' + target + '_' + template
        if not os.path.exists(directory):
            os.makedirs(directory)

        console('cp ./modeller_process.py ' + directory + '/modeller_process.py')

        process = Popen(["python3", "./modeller_process.py",
                         target, template, cwd], stdout=PIPE, stderr=PIPE, cwd=directory)
        stdout, stderr = process.communicate()

        console("rm -r " + directory)

        if len(stderr) > 0 or process.returncode != 0:
            print("MODELLER ERROR:", stderr)
        else:
            print("new model:", target, template)
            return ("ok", time.time() - start)
    return ("bad", 0)

processes_number = 48
pipe = Pipeline(UnorderedStage(run_modeller, processes_number))

waiting_candidates = 0
processed_candidates = 0
sum_time = 0

def get_pipe():
    global pipe
    global waiting_candidates
    global sum_time
    global processed_candidates
    res = pipe.get()
    waiting_candidates -= 1
    if res[0] == "ok":
        processed_candidates += 1
        sum_time += res[1]
    print("QUEUE:", waiting_candidates,
        "time:", sum_time / max(1, processed_candidates) / processes_number)
    return

import sys
import traceback

def process_id(id):
    global excluded_models
    global pipe
    global waiting_candidates
    global sum_time
    global processed_candidates
    print("START PROCESSING TARGET: ", id)
    get_pdb(id, False)
    candidates = get_candidates(id, get_chains(id)).split("|")
    target = id + get_chains(id)[0]
    for candidate_str in candidates:
        candidate = candidate_str.split(" ")
        if len(candidate) < 2:
            continue
        template = candidate[0] + candidate[1]
        try:
            if target == template:
                continue
            if (target, template) in excluded_models:
                #print("CACHE:", target, template)
                continue
            get_pdb(target, False)
            get_pdb(template, False)
            #print("Put candidate:", target, template)
            pipe.put(target + " " + template)
            waiting_candidates += 1
            print("QUEUE:", waiting_candidates,
                  "time:", sum_time / max(1, processed_candidates) / processes_number)
            while threading.active_count() > 1000:
                time.sleep(1)
            t = Thread(target=get_pipe)
            t.start()

        except urllib.error.HTTPError as e:
            print("HTTP: ", e.url)
        except Exception as e:
            print("EXCEPTION: ", e)
    return

class BaseThread(Thread):
    def __init__(self, callback=None, callback_args=None, *args, **kwargs):
        target = kwargs.pop('target')
        super(BaseThread, self).__init__(target=self.target_with_callback, *args, **kwargs)
        self.callback = callback
        self.method = target
        self.callback_args = callback_args

    def target_with_callback(self, *args):
        self.method(*args)
        if self.callback is not None:
            self.callback(*self.callback_args)


from threading import Lock
import threading
import queue

class MTQ():
    def __init__(self, target, max_threads):
        self.mutex = Lock()
        self.target = target
        self.max_threads = max_threads
        self.all_threads = []
        self.threads = set()
        self.tasks = queue.Queue()

    def cb(self):
        self.mutex.acquire()
        current_thread = threading.current_thread()
        if current_thread in self.threads:
            self.threads.remove(current_thread)

        while len(self.threads) < self.max_threads and self.tasks.qsize() > 0:
            t = BaseThread(target=self.target,
                           args=(self.tasks.get(),),
                           callback=self.cb,
                           callback_args=())
            print("LOAD REMAINING:", self.tasks.qsize())
            self.threads.add(t)
            self.all_threads.append(t)
            t.start()

        self.mutex.release()

    def add(self, arg):
        self.mutex.acquire()
        self.tasks.put(arg)
        self.mutex.release()
        self.cb()

    def join(self):
        for t in self.all_threads:
            t.join()

    def kill(self):
        for t in self.all_threads:
            t.terminate()


def main():
    global excluded_models
    global processes_number
    clean()
    excluded_models = get_all_generated_models()

    print("All models:", len(excluded_models))

    time.sleep(5)

    targets = []
    with open('./data/sids.txt', 'r') as file:
        for id in file.read().split('\n'):
            targets.append(id)
    if processes_number == 12:
        targets.reverse()

    mtq = MTQ(process_id, processes_number)

    for target in targets:
        mtq.add(target)

    mtq.join()

    pipe.put(None)

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
