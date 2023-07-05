import csv
import numpy as np

# define constants

ACTIVITE_LENGTH = 3600



def fish_input_process(file,time):
    row_data = []
    with open(file, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            row_data.append(row)
    row_data = np.array(row_data).astype(int)

    # Change row_data to float

    # Get first 10 columns
    row_data = row_data[:, 0:time]

    # Change all 0 to -1
    row_data[row_data == 0] = -1
    return row_data


if __name__ == '__main__':
    # file = input('Enter the name of the file: ')
    file = 'samples.csv'
    data = fish_input_process(file)
    print(data)
