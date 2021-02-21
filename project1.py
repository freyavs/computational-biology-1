from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
from itertools import combinations

f = open("woordenboek.txt", "r")
dictionary = f.read()
dictionary = dictionary.split('\n')
f.close()

dictionary = [ word.upper() for word in dictionary]

words = set() 
window_size = 0
for word in dictionary:
    # als het woord nummers of spaties heeft, kan het zeker niet in het aminozuur voorkomen
    if re.match('^[a-zA-Z]+$', word):
        words.add(word.upper())
        window_size = max(window_size, len(word))

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

    r = get_all_substring_combinations(s, res)
    return r

def get_all_substring_combinations(s, l):
    res = set()
    for word in l:
        # get all substrings that are not length 1
        substrings = (word[x:y] for x, y in combinations( 
            range(len(word) + 1), r = 2) if 1 < abs(x-y) < window_size)
        for w in substrings:
            if w in words: res.add(w)
            if w[::-1] in words: res.add(w[::-1])  
    
    return get_all_letters(s, res)

def get_all_letters(w,l):
    return l.union(set(w).intersection(words))


data = SeqIO.parse('covid.fasta', 'fasta')
for record in data:
    print(len(get_all_combinations(record)))

exit()

combs = get_all_combinations("TESTBESTAND")
print(combs)
print(len(combs))
exit()

#erge rna gif