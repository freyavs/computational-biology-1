data = SeqIO.parse('data/covid.fasta', 'fasta')
print(data)
for i, record in enumerate(data):
    f = open(f"outputs/eiwit_{i}.txt", "w+")
    combs = list(get_all_combinations(record))
    for c in sorted(combs):
        f.write(f'{c}\n')
    f.close()
    
exit()
    def search(self, key): 
        pCrawl = self.root 
        length = len(key) 
        for level in range(length): 
            letter = key[level]
            if not pCrawl.children.get(letter, None): 
                return False
            pCrawl = pCrawl.children[letter] 
  
        return pCrawl != None and pCrawl.string 