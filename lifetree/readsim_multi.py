#!/usr/bin/env python3

import sys
import os
import shutil
import subprocess
import random
from Bio.Seq import reverse_complement
from multiprocessing import Pool, set_start_method

readsim = '/public/lihaosen/anaconda3/envs/python2.7/bin/python /public/lihaosen/software/readsim-1.6/src/readsim.py'
epath = os.environ['PATH']

def runcmd(cmd, log, env=None):
    if env:
        os.environ['PATH'] = env + ':' + epath
        log.write('PATH: +' + env + "\n")
    log.write(' '.join(cmd) + "\n")
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while p.poll() is None:
        cmdout = p.stdout.readline().decode("utf8")
        if cmdout:
            log.write(cmdout)

    # restore to the original environment
    os.environ['PATH'] = epath

    # if error occurs, stop the script and exit
    if p.returncode != 0:
        log.write("\nError in '" + ' '.join(cmd) + "'!\nDirectly use the sequence as read!\n")
        return False
    else:
        return True

def run_mp(function, args_list, cpus, kwds={}):
    total = len(args_list)
    results = []

    #try:
    #    set_start_method('spawn')
    #except RuntimeError:
    #    pass
    p = Pool(cpus)
    finish = 0
    for args in args_list:
        if kwds:
            results.append(p.apply_async(function, args=args, kwds=kwds))
        else:
            results.append(p.apply_async(function, args=args))
    while True:
        finish_task = sum(1 for result in results if result.ready())
        if finish_task == finish:
            continue
        finish = finish_task
        sys.stdout.write("[{}] {}{}/{} ({}%)\r".format(("+" * int((finish/total)*50)) + (" " * (50 - int((finish/total)*50))), (" " * (len(str(total))-len(str(finish)))), finish, total, '%.2f' %(finish/total*100)))
        sys.stdout.flush()
        if finish == total:
            break
    p.close()
    p.join()
    results = [result.get() for result in results]

    # end the pbar
    print("\n")
    return results

def mp_harvest(dir_name, suffix, output):
    outfile = open(output, 'w')
    for file in os.listdir(dir_name):
        if file.endswith(suffix):
            for line in open(os.path.join(dir_name, file)):
                outfile.write(line)
    outfile.close()

def run_readsim(genome, tech, cov, output):
    log = open(output + '.log', 'w')
    cmd = readsim.split(' ')
    cmd.extend(['sim', 'fa', '--ref', genome, '--pre', output, '--tech', tech, '--cov_mu', cov, '--rev_strd', 'on'])
    ifcomplish = runcmd(cmd, log)
    if not ifcomplish:
        seqstr = ''
        seqid = None
        for line in open(genome):
            if line.startswith('>'):
                seqid = line.rstrip().lstrip('>').split(' ')[0]
            elif seqid:
                seqstr += line.rstrip()
        outfile = open(output + '.fasta', 'w')
        Fnum = 0
        Rnum = 0
        for i in range(int(cov)):
            if random.randint(0, 1) == 0:
                outfile.write(">{}:{}:F\n{}\n".format(seqid, i+1, seqstr))
                Fnum += 1
            else:
                outfile.write(">{}:{}:R\n{}\n".format(seqid, i+1, reverse_complement(seqstr)))
                Rnum += 1
        outfile.close()
        log.write("Directly generating reads: {} forward, {} reversed\n".format(Fnum, Rnum))
    log.close()

cpu = 8
if len(sys.argv) == 1 or sys.argv[1] == '-h':
	print("Usage: {} genome(fasta) output_prefix tech coverage cpu(default=8)".format(sys.argv[0]))
	sys.exit(0)
elif len(sys.argv) < 5:
	print("Error: options < 4!\nUsage: {} genome(fasta) output_prefix tech coverage cpu(default=8)".format(sys.argv[0]))
	sys.exit(1)
elif len(sys.argv) > 5:
	cpu = int(sys.argv[5])

if os.path.exists(sys.argv[2] + '.readsim.temp'):
	shutil.rmtree(sys.argv[2] + '.readsim.temp')
os.mkdir(sys.argv[2] + '.readsim.temp')

seqids = []
seqid = None
for line in open(sys.argv[1]):
	if line.startswith('>'):
		seqid = line.rstrip().lstrip('>').split(' ')[0]
		seqids.append(seqid)
		outfile = open(os.path.join(sys.argv[2] + '.readsim.temp', seqid + '.fna'), 'w')
	if seqid:
		outfile.write(line)	
outfile.close()

args_list = []
seqids.sort(key = lambda x : os.path.getsize(os.path.join(sys.argv[2] + '.readsim.temp', x + '.fna')), reverse = True)
for seqid in seqids:
	args_list.append((os.path.join(sys.argv[2] + '.readsim.temp', seqid + '.fna'), sys.argv[3], sys.argv[4], os.path.join(sys.argv[2] + '.readsim.temp', seqid)))
print("Running readsim in multiprocess...")
run_mp(run_readsim, args_list, cpu)
mp_harvest(sys.argv[2] + '.readsim.temp', '.fasta', sys.argv[2] + '.fasta')
mp_harvest(sys.argv[2] + '.readsim.temp', '.log', sys.argv[2] + '.readsim.log')
shutil.rmtree(sys.argv[2] + '.readsim.temp')

