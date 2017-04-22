import os

mainpath = os.path.dirname(__file__)
scratchpath = os.path.join(mainpath, "scratch")

def clear_ck(name):
    os.system("rm scratch/{}*.ck", name)

def clear_xml(name):
    os.system("rm scratch/{}*.xml", name)
