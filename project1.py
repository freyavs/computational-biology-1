from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from itertools import combinations

f = open("woordenboek.txt", "r")
words = f.read()
words = words.split('\n')
f.close()

words = set([ word.upper() for word in words])

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
            range(len(word) + 1), r = 2) if not abs(x-y) == 1)
        for w in substrings:
            if w in words or w[::-1] in words: res.add(w)
    
    return get_all_letters(s, res)

def get_all_letters(w,l):
    return l.union(set(w).intersection(words))


def get_all_reverse_combinations(l):
    l.update([word[::-1] for word in l])
    return l


def get_all_matches(s1, s2):
    return s1.intersection(s2)

data = SeqIO.parse('covid.fasta', 'fasta')
for record in data:
    print(get_all_combinations(record))

exit()

combs = get_all_combinations("TESTBESTAND")
print(combs)
exit()


#erge rna gif