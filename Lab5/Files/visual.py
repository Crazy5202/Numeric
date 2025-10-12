import matplotlib.pyplot as plt
import os

DATA_FOLDER = "results"

PATH = os.path.join(os.path.split(os.path.realpath(__file__))[0], DATA_FOLDER)

data_files = [file for file in os.listdir(PATH) if len(file)>1 and file[0]=='t' and file[1].isdigit() and file.endswith(".txt")]

for file in [data_files[0], data_files[len(data_files)//2], data_files[-1]]:
    full_path = os.path.join(PATH, file)

    x = []; u_solved = []; u_true = []
    with open(full_path, "r") as f:
        lines = f.readlines()
        for line in lines:
            line_stripped = line.strip().split(' ')
            x.append(float(line_stripped[0]))
            u_solved.append(float(line_stripped[1]))
            u_true.append(float(line_stripped[2]))
        plt.figure()
        plt.title("График в момент времени t=" + file.removesuffix('.txt')[1:])

        plt.plot(x, u_solved, marker='o', linestyle='-', color='red', label='solved')
        plt.plot(x, u_true, marker='o', linestyle='--', color='blue', label='true')

        plt.xlabel('x')
        plt.ylabel('u')
        
        plt.legend()
        plt.grid()
        plt.show(block=False)

while plt.get_fignums():
    plt.pause(1)