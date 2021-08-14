import os

def dataset_prep(fname = ''):
    if fname == '':
        print('please provide a file name or path!!!')
    else:
        dir = './datasets/'+fname
        my_file = open(dir)
        string_list = my_file.readlines()
        my_file.close()

        new_list=[]
        for line in string_list:
            nums = line.split()
            if len(nums)==1:
                new_list.append(nums[0]+'\n')
                continue
            nums[2] = ''
            new_list.append(' '.join(nums)+'\n')

        my_file = open(dir, "w")
        new_file_contents = "".join(new_list)

        my_file.write(new_file_contents)
        my_file.close()
    

paths = os.listdir('./datasets')
for dir in paths:
    print("preparing dataset in: ",dir)
    dataset_prep(dir)

