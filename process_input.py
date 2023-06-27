import csv
import numpy as np

#define constants

#ACTIVITE_LENGTH = 3600
ACTIVITE_LENGTH = 10




def fish_input_process(file):
    row_data = []
    with open(file, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            row_data.append(row)
    row_data = np.array(row_data)

    #Get first 10 columns
    row_data = row_data[:, 0:ACTIVITE_LENGTH].astype(int)

    #Change all 0 to -1
    row_data[row_data == 0] = -1
    return row_data




if __name__ == '__main__':
    #file = input('Enter the name of the file: ')
    file = 'fish1.csv'
    data = fish_input_process(file)