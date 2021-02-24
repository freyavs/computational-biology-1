from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
from dict_trie import *
import sys

def get_trie(dictionary_file):
    f = open(dictionary_file, "r")
    dictionary = ( word.upper() for word in f.read().split('\n') )
    f.close()
    words = Trie() 
    max_length = 0
    for word in dictionary:
        # als het woord nummers of spaties heeft, kan het zeker niet in het aminozuur voorkomen
        if re.match('^[A-Z]+$', word):
            words.add(word)
            max_length = max(max_length, len(word))

    return (words, max_length)

# handig om te testen met gewone strings
def fix_args(f):
    def aangepast(*args, **kwargs):
        new_args = [ SeqRecord(Seq(arg)) if isinstance(arg, str) else arg for arg in args ]
        return f(*new_args, **kwargs)
    return aangepast

@fix_args
def get_all_combinations(data):
    data = str(data.seq)
    res = set()
    max_jump = len(data) 
    res.add(data) #initieel woord met springlengte = 1
    count = 0
    for i in range(2,max_jump):
        for start_index in range(i):
            string = data[start_index::i]
            #sla deze op in een set zodat je al een groot deel duplicaten weghaald
            if len(string) > 1: 
                res.add(string)

    return res 

def get_all_substring_combinations(s, l, words, max_len):
    res = set()
    for word in l:
        length = len(word)
        for i in range(length):
            words.reset()
            dir_normal = words.find_update(word[i], 0)[1]
            dir_reverse = words.find_update(word[length-1-i], 1)[1]
            if dir_normal or dir_reverse:
                for j in range(i+2, length+1):
                    if j-i < max_len:
                        if dir_normal:
                            r = words.find_update(word[j-1], 0)
                            if r[0]: res.update(r[0])
                            dir_normal = r[1]

                        if dir_reverse:
                            r_reverse = words.find_update(word[length-j],1) 
                            if r_reverse[0]: res.update(r_reverse[0])
                            dir_reverse = r_reverse[1]

                        # als in geen van beide nog gezogd kan worden
                        if not dir_normal and not dir_reverse: break 

    # letters appart behandelen, die kunnen snel door in de speciale karakters te zoeken
    for letter in set(s):
        letters_found = [letter] + words.special_chars.get(letter, [])
        res.update(letters_found) 

    res.update(['X', 'O'])

    return res 

if __name__ == "__main__":
    trie, max_length = get_trie(sys.argv[1]) 

    data = SeqIO.parse(sys.argv[2], 'fasta')
    for record in data:
        print(f">{record.name}")
        combs = get_all_combinations(record)
        res = get_all_substring_combinations(record, combs, trie, max_length)
        res = sorted(list(res))
        for word in res:
            print(word)