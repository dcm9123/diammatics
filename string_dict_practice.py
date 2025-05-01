#Daniel Castaneda Mogollon, PhD
#Exercise 1 with dictionaries
'''
Exercise: Word Frequency Counter
Write a Python program that takes a string (a sentence or paragraph), and counts how many times each word appears. Store the results in a dictionary where the key is the word and the value is the count.

Instructions:
Ask the user to input a string.
Convert the string to lowercase (so The and the are considered the same).
Remove any punctuation (optional, for simplicity).
Split the string into words.
Count the frequency of each word using a dictionary.
Print the dictionary.
'''

input_string = input("Please provide a sentence of your choice:")
string_dict={}
new_list=[]
counter = 0
if isinstance(input_string, str):
    input_string = input_string.lower()
    words = input_string.split(" ")
    for word in words:
        for punctuation in '.,;!?"-_#$+=:':
            word = word.replace(punctuation,'')
        new_list.append(word)
    for items in new_list:
        if items in string_dict.keys():
            string_dict[items] +=1
        else:
            string_dict[items] = 1
    print("This is your modified sentence: "+" ".join(new_list))
else:
    raise TypeError("Please provide a sentence of your choice")
print(string_dict)