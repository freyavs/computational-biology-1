import numpy as np

class Node: 
      
    def __init__(self): 
        self.children = [None]*26
        self.string = None
  
class Trie: 
      
    def __init__(self): 
        self.special_chars = { 
            'N': [ 'B'],
            'D': ['B'],
            'Q': [ 'Z'],
            'E': [ 'Z'],
            'I': [ 'J'],    
            'V': ['U']
        }    
        self.root = Node() 
        self.search_nodes = [ [self.root], [self.root] ]

    def reset(self):
        self.search_nodes = [ [self.root], [self.root] ]
  
    def _charToIndex(self,ch): 
        return ord(ch)-ord('A') 
  
  
    def insert(self,key): 
        pCrawl = self.root 
        length = len(key) 
        for level in range(length): 
            index = self._charToIndex(key[level]) 
            if not pCrawl.children[index]: 
                pCrawl.children[index] = Node() 
            pCrawl = pCrawl.children[index] 
  
        pCrawl.string = key

    def search_update(self, letter, search_index):
        new_search_nodes = []
        results = []
        search = self.search_nodes[search_index]
        chars = [letter, 'X', 'O']  + self.special_chars.get(letter, []) 

        while self.search_nodes[search_index]:
            current_node = search.pop() 
            pset = list(filter(None, [ current_node.children[self._charToIndex(c)] for c in chars ]))
            new_search_nodes += pset 
            if pset: results += list(filter(None, [node.string for node in pset])) 

        if not new_search_nodes: return ([], False)
        self.search_nodes[search_index] = new_search_nodes
        return (results, True) 


    def search(self, key): 
        pCrawl = self.root 
        length = len(key) 
        for level in range(length): 
            index = self._charToIndex(key[level]) 
            if not pCrawl.children[index]: 
                return False
            pCrawl = pCrawl.children[index] 
  
        return pCrawl != None and pCrawl.isEndOfWord 
