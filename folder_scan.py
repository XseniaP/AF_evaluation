import os

def scan(directory):
    my_list = []
    for filename in os.listdir(directory):
        name = filename.split(".")[0]
        if "_" in name:
            my_list.append(name)
    print(my_list)
    return my_list

scan("/Users/kseniapolonsky/Downloads/RES")