import subprocess
import os

LIST = ["1ADQ_L","1CL7_H","1CL7_I","1CL7_L","1DEE_A","1DEE_B","1DFB_H", "1DFB_L", "1DN0_A","1DN0_B","1E6J_H","1E6J_L","1E6O_H","1E6O_L","1FN4_A","1FN4_B","1IL1_A","1IL1_B","1MVU_A","1MVU_B","1NGQ_H","1NGQ_L","1RHH_A","1RHH_B","1RZ8_A","1RZ8_B","1RZG_A","1RZG_B","1ZTX_H","1ZTX_L"]
# LIST = ["4WHT","4YNY","4ZTP","5E8E","5F6I","5FCU","5IHU","5IT2","5KW9","5M63","5XMH","6E9H","6OE5","6UE9","6VRY","7BWJ","7BXV","7CHF","7K75","7LXZ","8FAB"]
# LIST = ["2AGJ_H", "2AGJ_L", "2AJ3_A", "2AJ3_B", "2B1H_H", "2B1H_L", "2BX5_A", "2FB4_H", "2FB4_L", "2GHW_B", "2HRP_H",
#         "2HRP_L", "2J6E_H", "2J6E_L", "2PR4_H", "2PR4_L", "2QHR_H", "2QHR_L", "2X7L_H", "2X7L_L", "3DGG_A", "3DGG_B",
#         "3LMJ_H", "3LMJ_L", "3M8O_H", "3M8O_L", "3MBX_H", "3MBX_L", "3MLY_H", "3MLY_L", "3MNV_A", "3MNV_B", "3NTC_H",
#         "3NTC_L", "3Q6F_A", "3Q6F_B", "3Q6G_H", "3Q6G_L", "3QEG_H", "3QEG_L", "4HBC_H", "4HBC_L", "4JO4_H", "4JO4_L",
#         "4K3D_H", "4K3D_L", "4Q2Z_H", "4Q2Z_L"]


def replace(chain, fasta_name):
    file = open("script_main.txt", "r")
    replacement = ""
    # using the for loop
    # for line in file:
    #     line = line.strip()
    #     changes = line.replace("$1", chain)
    #     changes = changes.replace("$2", fasta_name)
    #     replacement = replacement + changes + "\n"

    data = file.read()
    # Searching and replacing the text
    # using the replace() function
    # path1 = "~/AP2/"+chain
    cwd = os.getcwd()
    # path1 = os.path.join(cwd, chain)
    path1 = "/groups/pupko/kseniap/AP2/" + chain
    # path2 = "~/AP2/"+fasta_name
    path2 = "/groups/pupko/kseniap/AP2/" + fasta_name
    data = data.replace("$1", path1)
    data = data.replace("$2", path2)
    file.close()

    # opening the file in write mode
    filename_new = "script_main2.txt"
    fout = open(filename_new, "w")
    fout.write(data)
    # fout.write(replacement)
    fout.close()


def run_AF():
    for chain in LIST:
        cwd = os.getcwd()
        fasta_name = "%s.fasta" % (chain)
        code = chain.split("_")[0]
        # create the folder
        path = os.path.join(cwd, chain)
        os.mkdir(path)
        replace(chain, fasta_name)

        rc = subprocess.Popen(["qsub", "script_main2.txt"])


if __name__ == '__main__':
    run_AF()
