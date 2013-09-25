#!/usr/bin/env python 

## \file autotest.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author Aniket C. Aranake, Alejandro Campos, Thomas D. Economon
#  \version 2.0.2.

import sys,time, os, subprocess, datetime, signal, os.path

class testcase:

  def __init__(self,tag_in):

    datestamp = time.strftime("%Y%m%d", time.gmtime())
    self.tag  = "%s_%s"%(tag_in,datestamp)  # Input, string tag that identifies this run

    # The test condition. These must be set after initialization
    self.test_iter = 1
    self.test_vals = []  

    # These can be optionally varied 
    self.case_dir    = "/home/ale11"
    self.joe_exec    = "joe"
    self.cfg_file    = "Joe.in"
    self.timeout     = 300
    self.tol         = 0.001
    self.outputdir   = "/home/ale11"

  def run_test(self):

    passed       = True
    exceed_tol   = False
    timed_out    = False
    iter_missing = True
    start_solver = True

    # Adjust the number of iterations in the config file   
    self.do_adjust_iter()

    # Assemble the shell command to run SU2
    self.joe_exec = os.path.join(os.environ['MUM_HOME'], 'bin', self.joe_exec)
    self.cfg_file = os.path.join('../inputfiles', self.cfg_file)
    command_base = "mpirun -np 8 %s %s > outputfile"%(self.joe_exec, self.cfg_file)
    command      = "%s"%(command_base)

    # Run Joe
    os.chdir(os.path.join(os.environ['MUM_HOME'], self.case_dir, 'bin')) 
    start   = datetime.datetime.now()
    process = subprocess.Popen(command, shell=True)  # This line launches Joe

    while process.poll() is None:
      time.sleep(0.1)
      now = datetime.datetime.now()
      if (now - start).seconds> self.timeout:
        try:
          process.kill()
          os.system('killall %s' % self.joe_exec)   # In case of parallel execution
        except AttributeError: # popen.kill apparently fails on some versions of subprocess... the killall command should take care of things!
          pass
        timed_out = True
        passed    = False

    # Examine the output
    f = open('outputfile','r')
    output = f.readlines()
    delta_vals = []
    sim_vals = []
    if not timed_out:
      start_solver = False
      for line in output:
        if not start_solver: # Don't bother parsing anything before Tecplot file is written
          if line.find('writeFlaggedCvsTecplot: full.000000.plt') > -1:
            start_solver=True
        else:   # Found the --Begin solver --- line; parse the input
          raw_data = line.split()
          try:
            iter_number = int(raw_data[1])
            data        = raw_data[6:]    # Take the last few columns for comparison
          except ValueError:
            continue
          except IndexError:
            continue
         
          if iter_number == self.test_iter:  # Found the iteration number we're checking for
            iter_missing = False
            if not len(self.test_vals)==len(data):   # something went wrong... probably bad input
              print "Error in test_vals!"
              passed = False
              break
            for j in range(len(data)):
              sim_vals.append( float(data[j]) )
              delta_vals.append( abs(float(data[j])-self.test_vals[j]) )
              if delta_vals[j] > self.tol:
                exceed_tol = True
                passed     = False
            break
          else:
            iter_missing = True

      if not start_solver:
        passed = False
        
      if iter_missing:
        passed = False

    print '=========================================================\n'
      
    # Write the test results 
    # for j in output:
    #  print j

    if passed:
      print "%s: PASSED"%self.tag
    else:
      print "%s: FAILED"%self.tag

    print 'execution command: %s'%command

    if timed_out:
      print 'ERROR: Execution timed out. timeout=%d'%self.timeout

    if exceed_tol:
      print 'ERROR: Difference between computed input and test_vals exceeded tolerance. TOL=%f'%self.tol

    if not start_solver:
      print 'ERROR: The code was not able to get to the "writeFlaggedCvsTecplot" section.'

    if iter_missing:
      print 'ERROR: The iteration number %d could not be found.'%self.test_iter

    print 'test_iter=%d, test_vals: '%self.test_iter,
    for j in self.test_vals:
      print '%f '%j,
    print '\n',

    print 'sim_vals: ',
    for j in sim_vals:
      print '%f '%j,
    print '\n',
  
    print 'delta_vals: ',
    for j in delta_vals:
      print '%f '%j,
    print '\n'
    
    return passed

  def do_adjust_iter(self):
  
    # Read the cfg file
    cfg_dir = os.path.join(os.environ['MUM_HOME'], self.case_dir, 'inputfiles')
    os.chdir(cfg_dir)
    file_in = open(self.cfg_file, 'r')
    lines   = file_in.readlines()
    file_in.close()
  
    # Rewrite the file with a .autotest extension
    self.cfg_file = "%s.autotest"%self.cfg_file
    file_out = open(self.cfg_file,'w')
    file_out.write('# This file automatically generated by autotest.py\n')
    file_out.write('# Number of iterations changed to %d\n'%(self.test_iter+1))
    for line in lines:
      if line.find("NSTEPS")==-1:
        file_out.write(line)
      else:
        file_out.write("NSTEPS = %d\n"%(self.test_iter+1))
    file_out.close()
    

if __name__=="__main__":
  '''This program runs SU^2 and ensures that the output matches specified values. This will be used to do nightly checks to make sure nothing is broken. '''

  # Build Joe
  os.system('make clean')
  os.system('make joe')

  os.chdir('./bin')
  if not os.path.exists("./joe"):
    print 'Could not build joe'
    sys.exit(1)
    
  os.chdir('../')

  # Flat Plate
  plateSA           = testcase('plateSA')
  plateSA.case_dir  = "vnv/flatplate"
  plateSA.cfg_file  = "fplate_sa.in"
  plateSA.test_iter = 100
  plateSA.test_vals = [3.1246e+03, 9.0161e-05]
  plateSA.joe_exec  = "joe 1"
  plateSA.timeout   = 1600
  plateSA.tol       = 0.001
  passed1            = plateSA.run_test()
  
  plateSST           = testcase('plateSST')
  plateSST.case_dir  = "vnv/flatplate"
  plateSST.cfg_file  = "fplate_sst.in"
  plateSST.test_iter = 100
  plateSST.test_vals = [3.0887e+03,4.8567e-03,7.4949e+06]
  plateSST.joe_exec  = "joe 2"
  plateSST.timeout   = 1600
  plateSST.tol       = 0.001
  passed2            = plateSST.run_test()
  
  plateEASMkom           = testcase('plateEASMkom')
  plateEASMkom.case_dir  = "vnv/flatplate"
  plateEASMkom.cfg_file  = "fplate_easmkom.in"
  plateEASMkom.test_iter = 100
  plateEASMkom.test_vals = [5.2915e+04, 3.0080e-02, 4.4128e+06]
  plateEASMkom.joe_exec  = "joe 6"
  plateEASMkom.timeout   = 1600
  plateEASMkom.tol       = 0.001
  passed3                = plateEASMkom.run_test()
  
  if (passed1 and passed2 and passed3):
    sys.exit(0)
  else:
    sys.exit(1)

