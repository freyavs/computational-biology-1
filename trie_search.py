from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
from itertools import combinations
from dict_trie import *

f = open("data/woordenboek.txt", "r")
dictionary = f.read()
dictionary = dictionary.split('\n')
f.close()

wordss = set()
words = Trie() 
window_size = 0
for word in dictionary:
    # als het woord nummers of spaties heeft, kan het zeker niet in het aminozuur voorkomen
    if re.match('^[a-zA-Z]+$', word):
        words.insert(word.upper())
        wordss.add(word.upper())
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
            length = len(string)
            if length > 1 : 
                res.add(string)

    r = get_all_substring_combinations(s, res)
    return r

def get_all_substring_combinations(s, l):
    res = set()
    for word in l:
        length = len(word)
        for i,letter in enumerate(word):
            words.reset()
            dir_normal = words.search_update(letter, 0)[1]
            dir_reverse = words.search_update(word[length-1-i], 1)[1]
            if dir_normal or dir_reverse:
                for j in range(i+1, length + 1):
                    #todo: check if i+1 can remove check, check if length+1 is necessary
                    if 1 < abs(i-j) < window_size:

                        if dir_normal:
                            r = words.search_update(word[j-1], 0)
                            if r[0]: res.update(r[0])
                            dir_normal = r[1]

                        if dir_reverse:
                            r_reverse = words.search_update(word[length-j],1) 
                            if r_reverse[0]: res.update(r_reverse[0])
                            dir_reverse = r_reverse[1]

                        if not dir_normal and not dir_reverse: break 
                        
    
    #return res
    return get_all_letters(s, res)

def get_all_letters(w,l):
    return l.union(set(w).intersection(wordss))



data = SeqIO.parse('data/covid.fasta', 'fasta')
for i, record in enumerate(data):
    combs = get_all_combinations(record)
    print(len(combs))


f.close()
    

data = SeqIO.parse('data/covid.fasta', 'fasta')
print(data)
for i, record in enumerate(data):
    f = open(f"outputs/eiwit_{i}.txt", "w+")
    combs = list(get_all_combinations(record))
    for c in sorted(combs):
        f.write(f'{c}\n')
    f.close()
    
exit()

combs = get_all_combinations("TESTEIWIT")
print(combs)
combs = get_all_combinations("TQEQSQTQEQIQWQIQT")
print(combs)
exit()
# tests



#erge rna gif