#!/usr/bin/python3
#coding:utf8

import time
import subprocess
from subprocess import call
import os
import sys

# import matplotlib
# matplotlib.use('Agg')

from tkinter import Tk
from tkinter.filedialog import askopenfilename

if "WD" in os.environ:
  os.chdir(os.environ["WD"])
  path="/app"
else:
  path=os.getcwd()
  
def main ():
  start= time.time()

  if len(sys.argv)>=2:
    parameter_file=sys.argv[1]
  else:
    Tk().withdraw() 
    parameter_file=askopenfilename(title="Choose parameter file")  

    
  if len(parameter_file)==0:
    print("Missing parameter file.")
    sys.exit(0)

    
  os.chdir(os.path.dirname(parameter_file))
  
  choice=input('Use DeCoSTAR ? (y/n) :')
  
  if choice=='y':                           
    f= path + "/DeCoSTAR parameter.file="+parameter_file  
    subprocess.call(f, shell=True)
  
  print("Diagnostic in process ...")
  subprocess.call(['python3', path + '/Diagnostic.py', parameter_file])

  
#  print("\n#analysis done in " + str(time.time() - start) + " sec")


if __name__=="__main__":
  main()
