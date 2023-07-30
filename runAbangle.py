import pathlib
import subprocess

import pandas as pd
PATH = '/Users/kpolonsky/PycharmProjects/rmsd/rmsd/'

def runList_Abangle(mylist):
    # df = pd.read_excel(io=filename, sheet_name='Sheet1')
    # mylist = df['pdb'].tolist()
    # molec = ""
    # for m in mylist:
    #     molec += PATH + m + '.pdb' + " "
    #     molec += PATH + m + '_ranked_0.pdb' + " "

    # mylist =['7l77_HL']
    for m in mylist:
        molec = PATH + m + '.pdb' + " "
        # args = "/ABangle" + " -i " + molec + " " + "-outdir /Users/kpolonsky/PycharmProjects/rmsd/rmsd/ -store csv -name Abangle"
        args = "/ABangle" + " -i " + molec + " " + "-outdir /Users/kpolonsky/PycharmProjects/rmsd/rmsd/ -store stdout -name Abangle"

        cmd = "/usr/bin/python3 ." + args
        pr = open(PATH + "/Abanglelog.txt", 'a')
        pr.write(molec + "\n")
        try:
            p = subprocess.Popen(cmd, stdout=pr, shell=True)
        except:
            print("An exception occurred " + molec)
            pr.write("An exception occurred " + molec + "\n")
        pr.close()

        # molec = PATH + m + '_ranked_0.pdb' + " "
        # args = "/ABangle" + " -i " + molec + " " + "-outdir /Users/kpolonsky/PycharmProjects/rmsd/rmsd/ -store csv -name Abangle"
        #
        # cmd = "/usr/bin/python3 ." + args
        # pr = open(PATH + "/Abanglelog.txt", 'a')
        # try:
        #     p = subprocess.Popen(cmd, stdout=pr, shell=True)
        # except:
        #     print("An exception occurred " + molec)
        #     pr.write(molec + "\n")
        # pr.close()

    # molec = '/Users/kpolonsky/PycharmProjects/rmsd/rmsd/7jlk_AB.pdb '
    # args = "/ABangle" + " -i " + molec + " " + "-outdir /Users/kpolonsky/PycharmProjects/rmsd/rmsd/ -store csv -name Abangle"

    # cmd = "/usr/bin/python3 ." + args
    # try:
    #     pr = open(PATH + "/Abanglelog.txt", 'a')
    #     p = subprocess.Popen(cmd, stdout=pr, shell=True)
    # except:
    #     print("An exception occurred " + molec)
    #     pr.write(molec + "\n")
    #     pr.close()

if __name__ == '__main__':
    # runList_Abangle('/Users/kpolonsky/PycharmProjects/rmsd/rmsd/HL_list.xlsx')
    runList_Abangle(['6WAS_HL_ImmuneBuilder_model'])