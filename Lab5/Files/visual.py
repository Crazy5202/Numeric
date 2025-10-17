import matplotlib.pyplot as plt
import os
from natsort import natsorted

# DATA_FOLDER = "results"
# PATH = os.path.join(os.path.split(os.path.realpath(__file__))[0], DATA_FOLDER)

def visualise(path):

    data_files = [file for file in os.listdir(path) if file.endswith(".txt") and file[0]!='p']

    data_files = natsorted(data_files)

    figure = plt.figure(figsize=(14, 9))
    counter = 1

    chosen_files = [data_files[len(data_files)//4], data_files[len(data_files)//2], data_files[-1]]

    for file in chosen_files:
        full_path = os.path.join(path, file)

        x = []; u_solved = []; u_true = []
        with open(full_path, "r") as f:
            lines = f.readlines()
            t_str = lines[0]
            for line in lines[1:]:
                line_stripped = line.strip().split(' ')
                x.append(float(line_stripped[0]))
                u_solved.append(float(line_stripped[1]))
                u_true.append(float(line_stripped[2]))
            
            p = figure.add_subplot(2,len(chosen_files),counter)

            p.plot(x, u_solved, marker='o', linestyle='-', color='red', label='solved')
            p.plot(x, u_true, marker='o', linestyle='--', color='blue', label='true')

            plt.title("График в момент времени t=" + t_str)
            plt.xlabel('x')
            plt.ylabel('u')
            plt.grid()
            plt.legend()
        counter += 1

    pogr_path = path + '/p.txt'
    with open(pogr_path, "r") as f:
        lines = f.readlines()
        x = []; pogr = []
        for line in lines:
            line_stripped = line.strip().split(' ')
            x.append(float(line_stripped[0]))
            pogr.append(float(line_stripped[1]))

        p = figure.add_subplot(2,1,2)

        p.plot(x, pogr, marker='', linestyle='-', color='red')

        plt.title("Погрешность в зависимости от времени")
        plt.xlabel('x')
        plt.ylabel('norm')
        plt.grid()

    plt.tight_layout()
    plt.show(block=True)

# visualise(PATH)