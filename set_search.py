from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
from itertools import combinations
import sys
from dict_trie import *


#dit geeft niet de gevraagde output want dit was gewoon een probeersel

def get_set(dict_file):
    f = open(dict_file, "r")
    dictionary = ( word.upper() for word in f.read().split('\n') )
    f.close()
    words = set() 
    max_len = 0
    for word in dictionary:
        # als het woord nummers of spaties heeft, kan het zeker niet in het aminozuur voorkomen
        if re.match('^[A-Z]+$', word):
            words.add(word)
            max_len = max(max_len, len(word))

    return words, max_len

def fix_args(f):
    def aangepast(*args, **kwargs):
        new_args = [ SeqRecord(Seq(arg)) if isinstance(arg, str) else arg for arg in args ]
        return f(*new_args, **kwargs)
    return aangepast

@fix_args
def get_all_combinations(s):
    s = str(s.seq)
    res = set()
    max_jump = len(s) 
    res.add(s) #initieel woord met springlengte = 1
    for i in range(2,max_jump):
        for start_index in range(i):
            string = s[start_index::i]
            # voeg enkel strings van lengte groter dan 1 toe
            if len(string) > 1: 
                res.add(string)

    return res

def get_all_substring_combinations(s, l, words, max_len):
    res = set()
    for word in l:
        # get all substrings that are not length 1
        substrings = (word[x:y] for x, y in combinations( 
            range(len(word) + 1), r = 2) if 1 < abs(x-y) < max_len)
        for w in substrings:
            if w in words: res.add(w)
            if w[::-1] in words: res.add(w[::-1])  

    # voeg letters appart toe 
    res.update(set(s).intersection(words))
    return res 

if __name__ == "__main__":
    words, max_length = get_set(sys.argv[1]) 

    data = SeqIO.parse(sys.argv[2], 'fasta')
    for record in data:
        print(record.name)
        combs = get_all_combinations(record)
        res = get_all_substring_combinations(record, combs, words, max_length)
        print(len(res))
