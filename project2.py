from trie2 import * 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
from itertools import combinations


def make_trie():
    f = open("woordenboek.txt", "r")
    dictionary = f.read()
    dictionary = dictionary.split('\n')
    f.close()

    words = Trie()
    window_size = 0
    for word in dictionary:
        # als het woord nummers of spaties heeft, kan het zeker niet in het aminozuur voorkomen
        if re.match('^[a-zA-Z]+$', word):
            words.insert(word.upper())
            window_size = max(window_size, len(word))

    return words, window_size

def fix_args(f):
    def aangepast(*args, **kwargs):
        new_args = [ SeqRecord(Seq(arg)) if isinstance(arg, str) else arg for arg in args ]
        return f(*new_args, **kwargs)
    return aangepast

@fix_args
def get_all_combinations(s, words, window_size):
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

    r = get_all_substring_combinations(s, res, words, window_size)
    return r

def get_all_substring_combinations(s, l, words, window_size):
    res = set()
    for word in l:
        # get all substrings that are not length 1
        substrings = (word[x:y] for x, y in combinations( 
            range(len(word) + 1), r = 2) if 1 < abs(x-y) < window_size)
        for w in substrings:
            if words.search(w): 
                res.add(w)
            if words.search(w[::-1]): 
                res.add(w[::-1])
    
    for letter in s:
        if words.search(w): res.add(w)
    
    return res


if __name__ == "__main__":
    words = make_trie()
    print("Trie built.")

    data = SeqIO.parse('covid.fasta', 'fasta')
    for record in data:
        print(" --------------- eiwit ------------------ ")
        print(get_all_combinations(record, words[0],words[1]))
    exit()
    print(words.search("TEST"))

    combs = get_all_combinations("TESTBESTAND", words)
    print(combs)
    exit()

    data = SeqIO.parse('covid.fasta', 'fasta')
    for record in data:
        print(" --------------- eiwit ------------------ ")
        print(get_all_combinations(record, words))
    exit()




#erge rna gif