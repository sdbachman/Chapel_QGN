import numpy as np
import subprocess
import re
import os


root_dir = os.getcwd()
print('Current directory is' + root_dir)

A_start = 2
A_end = 11
A_int = 1
pows = np.linspace(A_start,A_end, A_end-A_start+1)
pows = np.arange(A_start, A_end+1, A_int)
A = 3*(2**pows)

B_start = 2
B_end = 20
B_int = 2
B = np.linspace(B_start,B_end, int( (B_end-B_start+2)/2) )
B = np.arange(B_start, B_end+1, B_int)

print("There are ", len(A), " horizontal domain sizes:")
print(A)
print()

print("There are ", len(B), " layers:")
print(B)
print()

for a in A:
  for b in B:

    a_str = str(int(a))
    b_str = str(int(b))

    # Make directory
    dirstr = a_str + 'x' + a_str + 'x' + b_str

    command = 'mkdir ' + dirstr
    mkd = subprocess.Popen(command, shell = True)
    mkd.wait()

    # Copy code
    command = 'cp ./code/* ./' + dirstr
    mkd = subprocess.Popen(command, shell = True)
    mkd.wait()

    # Enter directory
    print('Changing to directory ' + dirstr)
    os.chdir(dirstr)

    # Set up strings to replace
    f = open('test.f90')
    text = f.read()
    f.close()
    SF = re.compile('REPL1')
    SF2 = re.compile('REPL2')
    SF3 = re.compile('REPL3')
    new = SF.sub('integer, parameter :: nx = ' + a_str,text)
    new2 = SF2.sub('integer, parameter :: ny = ' + a_str,new)
    new3 = SF3.sub('integer, parameter :: nz = ' + b_str,new2)
    f = open('test.f90', 'w')
    f.write(new3)
    f.close()

    # Submit job
    command = 'qsub run_test.sh'
    mkd = subprocess.Popen(command, shell = True)
    mkd.wait()

    print(" ")
    os.chdir(root_dir)

