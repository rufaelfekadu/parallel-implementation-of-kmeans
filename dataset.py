import pandas as pd

my_file = open("dataset-10000.txt")
string_list = my_file.readlines()
my_file.close()
new_list=[]
for line in string_list:
    nums = line.split()
    new_list.append(','.join(nums)+'\n')

new_list[0] = "x,y,z\n"

my_file = open("data.txt", "w")
new_file_contents = "".join(new_list)

my_file.write(new_file_contents)
my_file.close()

data = pd.read_csv('data.txt')
dataset = data.drop('z',axis = 1)
dataset.to_csv('dataset.csv',index=False)
