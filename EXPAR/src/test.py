import todo
import elixir
import configure_finger
import TmDeltaG
import SecStructures_jf5
from Bio.SeqRecord import SeqRecord
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import subprocess
import os
def main():
    try:
        a=os.popen('del sowien.txt >junk.txt').read()
        os.popen('del sowien.txt')
    except:
        print "ok"
    print a
if __name__ == "__main__":
    main()