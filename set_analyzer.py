import argparse

parser = argparse.ArgumentParser(description="This program was created to determine the items that are shared between two sets of files with only the"
                                 "item in each line")
print('''This job requires only the path to each individual txt file and the name of your output file:
--file1 [STR] [Path to the first txt file]
--file2 [STR] [Path to the second txt file]
--out_file [STR] [Name and path of your output file]
''')

parser.add_argument('--file1', type=str)
parser.add_argument('--file2', type=str)
parser.add_argument('--out_file', type=str)

args = parser.parse_args()
file1 = args.file1
file2 = args.file2
out_file = args.out_file

list1=[]
list2=[]
output_file = open(out_file, "w")
with open(file1, "r") as file:
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        list1.append(line)
with open(file2, "r") as file:
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        list2.append(line)
        
unique_list1 = list(set(list1) - set(list2))
unique_list2 = list(set(list2) - set(list1))
intersection = list(set(list1) & set(list2))
union = list(set(list1) | set(list2))

print("There are "+str(len(unique_list1))+" items present uniquely in your file 1 "+str(file1)+" are:")
for items in unique_list1:
     print(items, end=',')
print("")
print("There are "+str(len(unique_list2))+" items present uniquely in your file 2 "+str(file2)+" are:")
for items in unique_list2:
     print(items, end=',')
print("")
print("The are "+str(len(intersection))+" items present in the intersection of both files:")
for items in intersection:
     print(items, end=',')
print("")
print("The are "+str(len(union))+" items present in the union of both files are:")
for items in union:
     print(items,end=',')
print("")

output_file.write("There are "+str(len(unique_list1))+" items present uniquely in your file 1 "+str(file1)+" are:")
for items in unique_list1:
     output_file.write(items)
     output_file.write(",")
output_file.write("\nThere are "+str(len(unique_list2))+" items present uniquely in your file 2 "+str(file2)+" are:")
for items in unique_list2:
     output_file.write(items)
     output_file.write(",")
output_file.write("\nThe are "+str(len(intersection))+" items present in the intersection of both files:")
for items in intersection:
     output_file.write(items)
     output_file.write(",")
output_file.write("\nThe are "+str(len(union))+" items present in the union of both files are:")
for items in union:
     output_file.write(items)
     output_file.write(",")
