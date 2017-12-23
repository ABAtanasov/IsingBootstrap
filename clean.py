import os
import sys

mainpath = os.path.dirname(__file__)
scratchpath = os.path.join(mainpath, "scratch")

def clear_ck(name=""):
    os.system("rm scratch/{}*.ck".format(name))

def clear_xml(name=""):
    os.system("rm scratch/{}*.xml".format(name))

if __name__ == "__main__":
    if sys.argv[1] == "ck":
        clear_ck()
    elif sys.argv[1] == "xml":
        clear_xml()

