import os

path = "/Users/danielcm/Desktop/"
os.chdir(path)
string_dict = {}
new_list = []
stopwords = [
    "the", "a", "in", "of", "and", "to", "for", "on", "with", "at", "by",
    "an", "be", "this", "that", "it", "from", "as", "is", "are", "was", "but"
]
with open("dict_exercise_2.txt","r") as in_file:
    lines = in_file.readlines()
    for line in lines:
        for word in line.split(" "):
            word = word.lower()
            for i in ".,/';:+=!?":
                word = word.replace(i, "")
            new_list.append(word)

for item in new_list:
    if item not in stopwords:
        if item in string_dict:
            string_dict[item] = string_dict[item] + 1
        else:
            string_dict[item] = 1
    else:
        continue

sorted_list = sorted(string_dict.items(), key = lambda x: x[1], reverse=True)
count = 0

print("The top five words with most occurrences are: ")
for item in sorted_list:
    print(str(item[0])+" is found "+str(item[1])+" times.")
    count = count + 1
    if count > 5:
        break


