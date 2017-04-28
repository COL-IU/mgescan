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
        if not cmd_exists("hmmsearch"):
            print ("Prerequisite software: HMMER is not installed.")
            print ("Please install HMMER v3.1b1 and re-run this setup.")
            print ("Exiting without installation.")
            sys.exit()
        if not cmd_exists("transeq"):
            print ("Prerequisite software: EMBOSS is not installed.")
            print ("Please install EMBOSS v6.6.0 and re-run this setup.")
            print ("Exiting without installation.")
            sys.exit()
        if not cmd_exists("trf"):
            print ("Prerequisite software: Tandem Repeats Finder is not installed.")
            print ("Please install Tandem Repeats Finder v407b and re-run this setup.")
            print ("Exiting without installation.")
            sys.exit()
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
            sys.exit()

        os.system("cp -pr * "+install_path)

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

