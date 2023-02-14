import numpy as np
import subprocess
import re
import os
import glob
import sys

first = int(sys.argv[1])
last = int(sys.argv[2])

def source(script, update=1):
    pipe = subprocess.Popen(". %s; env" % script, stdout=subprocess.PIPE, shell=True)
    data = pipe.communicate()[0]
    env = dict((line.split("=", 1) for line in data.splitlines()))
    print(env)
    if update:
        os.environ.update(env)
    return env

root_dir = os.getcwd()
print('Current directory is' + root_dir)

base_dir = "/glade/scratch/bachman/Ian_QG/QGN_Beta/TEST/CHPL/FFT_SPEED_TEST/FORTRAN/"
dirs = glob.glob(base_dir + "[1-9]*")

dd = []
for d in dirs:
  dd.append( d.split("/")[-1] )

print("Making directories: ")
print(dd[first:last])

for d in dd[first:last]:
# Make directory

    sizes = d.split("x")
    nz = sizes[-1]
    nx = sizes[0]

    command = 'mkdir ' + d
    mkd = subprocess.Popen(command, shell = True)
    mkd.wait()

    # Copy code
    command = 'cp ./code/* ./' + d
    mkd = subprocess.Popen(command, shell = True)
    mkd.wait()

    # Enter directory
    print('Changing to directory ' + d)
    os.chdir(d)

    # Submit job
    #command = 'source ~/.chapel_pbs-gasnet'
    #mkd = subprocess.Popen(command, shell = True)
    #mkd.wait()


    # Set up strings to replace
    f = open('run.sh')
    text = f.read()
    f.close()
    SF = re.compile('REPL1')
    command = './test -nl 1 --filename=' + base_dir + d + "/test_grid" + " --nx=" + str(nx) + " --ny=" + str(nx) + " --nz=" + str(nz) # + " &"
    print(command)
    new = SF.sub(command,text)
    f = open('run.sh', 'w')
    f.write(new)
    f.close()

    env = {}
    env.update(os.environ)

    command = './run.sh'
    #command = './test -nl 1 --filename=' + base_dir + d + "/test_grid" + " --nx=" + str(nx) + " --ny=" + str(nx) + " --nz=" + str(nz)
    #mkd = subprocess.Popen(command, shell = True)
    #mkd = subprocess.Popen(command, shell = True, env=env)
    #mkd.wait()

    os.system(command)

    print(" ")
    os.chdir(root_dir)
