import sys
import subprocess
import shlex
import glob
import os
import errno
import shutil
import datetime
import re

def myLog(msg,f=sys.stdout):
 now = str(datetime.datetime.now()).split('.')[0]
 f.write("\n" + now + "\t" + msg + "\n")
 f.flush()
 
def cleanDirPath(dirPath):
 return os.path.abspath(dirPath) + "/"
 
def cleanFilePath(filePath):
 return os.path.abspath(filePath)
 
def makeDir(path):
 path = cleanDirPath(path)
 try:
  os.makedirs(path)
 except OSError as exception:
  if exception.errno != errno.EEXIST:
   raise

def getDirPath(filePath):
 return os.path.abspath(os.path.dirname(filePath)) + "/"
   
def getWorkingDir():
 return os.path.abspath(os.path.dirname(sys.argv[0])) + "/"
 
def makeTempDir(outdir):
 pid = os.getpid()
 tempdir = cleanDirPath(outdir) + str(pid) + "/"
 makeDir(tempdir)
 return tempdir
 
def removeDir(path):
 shutil.rmtree(path)
 
def runCommand(cmd,log=sys.stdout):
 if log != None:
  myLog("Running command "+cmd,log)
 p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 out, err = p.communicate()
 return out, err
 
def runCommandShell(cmd,log=sys.stdout):
 if log != None:
  myLog("Running command "+cmd,log)
 p = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
 out, err = p.communicate()
 return out, err
 
def listFiles(indir,wildcard='*'):
 files = glob.glob(os.path.join(indir, wildcard))
 fullnames = map(os.path.basename,files)
 return fullnames
 
def removeExt(filename):
 return os.path.splitext(filename)[0]
 
def myljust(val,offset):
 if type(val) is str:
  return val.ljust(offset)
 elif type(val) is float:
  return "{0:.2f}".format(val).ljust(offset)
 elif type(val) is int:
  return str(val).ljust(offset)

def find_between(s, start, end):
 result = re.search(start+'(.*?)'+end, s)
 if result == None:
  return None
 else:
  return result.group(1)
