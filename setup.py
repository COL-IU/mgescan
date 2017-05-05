import os
import sys
import time
import subprocess as sub
from setuptools import setup
from setuptools.command.bdist_egg import bdist_egg
def cmd_exists(cmd):
    return sub.call(["which", cmd], stdout=sub.PIPE, stderr=sub.PIPE) == 0

class MGEScanInstall(bdist_egg):
    def run(self):

        '''
        Checking Prerequisites
        '''
        hmmer_exists = cmd_exists("hmmsearch")
        emboss_exists = cmd_exists("transeq")
        trf_exists = cmd_exists("trf")
        blast_exists = cmd_exists("blastn")
        if not (hmmer_exists and emboss_exists and trf_exists and blast_exists):
            print ("[Error] Prerequisite software(s):")
            if not hmmer_exists:
                print ("\t- HMMER")
            if not emboss_exists:
                print ("\t- EMBOSS")
            if not trf_exists:
                print ("\t- Tandem Repeats Finder")
            if not blast_exists:
                print ("\t- NCBI-BLAST+")
            print ("\tare not installed!")
            print ("Please install the above software and re-run this setup.")
            print ("Exiting without installation!")
            sys.exit()

        '''
        Checking Install Path
        '''
        install_path = ""
        if not os.environ.get('MGESCAN_HOME'):
            print ("$MGESCAN_HOME is not defined where MGESCAN will be" + \
                    " installed.")
            def_home = raw_input("Would you like to install MGEScan at " + \
                    os.environ.get("HOME") + "/mgescan4 (Y/n)?")
            if def_home.lower() == 'n':
                print ("Run 'export MGESCAN_HOME=<your desired destination" + \
                        " path to install>' if you want to install somewhere"+\
                        " else\n")
                sys.exit()
            install_path = os.environ.get('HOME') + "/mgescan4"
        else:
            print ("Installing MGEScan in "+os.environ.get('MGESCAN_HOME')+" (MGESCAN_HOME)")
            install_path = os.environ.get('MGESCAN_HOME')
        if not os.path.exists(install_path):
            os.makedirs(install_path, 0755)
        else:
            print ("It looks like a previous installation of MGEScan already exists in "+install_path)
            print ("Please remove the previous installation and re-run this setup.")
            print ("Exiting without installation!")
            sys.exit()

        '''
        Installing
        '''
        if cmd_exists("make") and cmd_exists('gcc') and cmd_exists('g++'):
            os.system("cd mgescan/ltr/MER; make clean; make")
            os.system("cd mgescan/nonltr/; make clean; make translate")
            os.system("cd mgescan/nonltr/hmm;make clean; make")
        else:
            print ("[Warning] make|gcc|g++ does not exist. Compile code is skipped")
            time.sleep(3)
        if cmd_exists("mpicc"):
            os.system("cd mgescan;mpicc mpi_mgescan.c -o mpi_mgescan")
        else:
            print ("[Warning] mpicc does not exist. Compile mpi code is"\
                    + " skipped")
            time.sleep(3)

        os.system("cp -pr * "+install_path)

        if not os.environ.get('GALAXY_HOME'):
            print ("[Warning] $GALAXY_HOME is not defined. Copying MGEScan Galaxy files is skipped.")
            time.sleep(3)
        else:
            os.system("cp -pr galaxy-modified/* $GALAXY_HOME")

        bdist_egg.run(self)

class MGEScanInstallOnly(bdist_egg):
    def run(self):
        bdist_egg.run(self)

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

reqs = [line.strip() for line in open('requirements.txt')]
setup(
        name = "MGEScan",
        version = "4.0.0",
        author = "Wazim Mohammed Ismail",
        author_email = "wazimoha@indiana.edu",
        description = ("MGEScan on Galaxy Workflow System for identifying transposable "
            "elements in genome sequences"),
        license = "GPLv3",
        keywords = "MGEScan, Galaxy workflow",
        url = "https://github.com/COL-IU/mgescan",
        packages = ['mgescan','mgescan.dna'],
        install_requires = reqs,
        long_description = read('README.md'),
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Scientific/Engineering",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU GEneral Public License v3 (GPLv3)",
            "Operating System :: POSIX :: Linux",
            "programming Language :: Python",
            ],
        entry_points='''
            [console_scripts]
            mgescan=mgescan.cmd:main
            mg_split=mgescan.split:main
            nonltr=mgescan.nonltr:main
            ''',

        cmdclass={'bdist_egg': MGEScanInstall},  # override bdist_egg
        )

