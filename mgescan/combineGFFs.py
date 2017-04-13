"""combineGFFs.py: combine GFFs from LTR, nonLTR and DNA

Usage:
	combineGFFs.py <datadir>
	combineGFFs.py (-h | --help)
	combineGFFs.py --version

Options:
	-h --help   Show this screen.
	--version   Show version.

"""
from docopt import docopt
import os
from mgescan import utils

class CombineGFFs(object):

	def __init__(self, args=None):
		self.args = args
		if args:
			self.set_parameters()

	def set_parameters(self):
		self.set_datadir(self.args['<datadir>'])

	def set_datadir(self, path):
		self.datadir = utils.get_abspath(path)
		return self.datadir

	def run(self):
		self.combine()

	def combine(self):
		datadir = self.datadir
		masterlist = []
		if os.path.isfile(datadir+"/ltr/ltr.gff3"):
			f = open(datadir+"/ltr/ltr.gff3")
			for line in f:
				masterlist.append(line.split(None))
			f.close()
		if os.path.isfile(datadir+"/info/nonltr.gff3"):
			f = open(datadir+"/info/nonltr.gff3")
			for line in f:
				masterlist.append(line.split(None))
			f.close()
		if os.path.isfile(datadir+"/dna/dna.gff3"):
			f = open(datadir+"/dna/dna.gff3")
			for line in f:
				masterlist.append(line.split(None))
			f.close()
		f = open(datadir+"/mge.gff3","w")
		sortedlist = sorted(masterlist,key=lambda x: (x[0],x[3]))
		for item in sortedlist:
			f.write("\t".join(item))
			f.write("\n")
		f.close()

def main():
	arguments = docopt(__doc__, version='combineGFFs.py 0.1')
	combineGFFs = CombineGFFs(arguments)
	combineGFFs.run()

if __name__ == "__main__":
	main()
