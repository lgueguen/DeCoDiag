#!/usr/bin/python3
#coding:utf8

import time
import subprocess
from subprocess import call
import os
import sys
import diagnostic

# import matplotlib
# matplotlib.use('Agg')

from tkinter import Tk
from tkinter.filedialog import askopenfilename

if "WD" in os.environ:
  os.chdir(os.environ["WD"])
  path="/app"
else:
  path=os.getcwd()
  


start= time.time()

if len(sys.argv)>=3:
  wd_dir=sys.argv[1]
  parameter_file=sys.argv[2]
elif len(sys.argv)==2:
  if os.path.isfile(sys.argv[1]):
    wd_dir=os.path.dirname(sys.argv[1])
    parameter_file=os.path.basename(sys.argv[1])
  else:
    try:
      Tk().withdraw() 
      parameter_file=askopenfilename(title="Choose parameter file")  
      wd_dir=os.path.dirname(parameter_file)
      parameter_file=os.path.basename(parameter_file)
    except:
      print("Window interface not available. Please entre a parameter file.")
      sys.exit(0)
    
if len(parameter_file)==0:
  print("Missing parameter file.")
  sys.exit(0)

if os.path.isdir(wd_dir):
    os.chdir(wd_dir)
else:
  os.chdir(os.path.dirname(wd_dir))

choice=input('Use DeCoSTAR ? (y/n) :')

if choice=='y':                           
  print(path + "/DeCoSTAR parameter.file="+parameter_file)  
  f= path + "/DeCoSTAR parameter.file="+parameter_file  
  subprocess.call(f, shell=True)


####################################################
#####

print("Diagnostic in process ...")

diag=diagnostic.Diagnostic(parameter_file)
diag.build_ancestral()

print(diag.anc)

diag.output_cycles()

  #   ######################
  # # zip 6-cycles with duplications

  # chzip=input('Zip trees with duplicated families? (y/n): ')

  # if chzip=='y':                   
  # new_directory=input('Name of the new directory for trees ? : ')
  # # while os.path.exists(new_directory):
  # #   print("This directory already exists.")
  # #   new_directory=input('Name of the new directory for trees ? : ')
  
  # zipfam=self.zip_cycles6_dup()

  # new_config_file=self.increment_suffix_in_param_file("gene.distribution.file")
  # self.output_gene_trees(zipfam, new_directory, new_config_file)
  # self.dfile["gene.distribution.file"]=new_config_file
  # self.new_param_file()

  
#  print("\n#analysis done in " + str(time.time() - start) + " sec")

