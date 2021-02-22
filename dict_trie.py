import numpy as np

class TrieNode: 
      
    # Trie node class 
    def __init__(self): 
        #self.children = [None]*26
        self.children = dict()
        self.string = None 
        # isEndOfWord is True if node represent the end of the word 
        self.isEndOfWord = False
        self.isleaf = False
  
class Trie: 
      
    # Trie data structure class 
    def __init__(self): 
        self.root = self.getNode() 
        self.search_nodes = [ [self.root], [self.root] ]

    def reset(self):
        self.search_nodes = [ [self.root], [self.root] ]
  
    def getNode(self): 
      
        # Returns new trie node (initialized to NULLs) 
        return TrieNode() 
  
    def insert(self,key): 
          
        # If not present, inserts key into trie 
        # If the key is prefix of trie node,  
        # just marks leaf node 
        pCrawl = self.root 
        length = len(key) 
        for level in range(length): 
            #index = self._charToIndex(key[level]) 
            index = key[level]
            # if current character is not present 
            if not pCrawl.children.get(index, None): 
                pCrawl.children[index] = self.getNode() 
            pCrawl = pCrawl.children[index] 
  
        # mark last node as leaf 
        pCrawl.string = key
        pCrawl.isEndOfWord = True

    def search_update(self, letter, search_index):
        special_chars = { 
            'N': ['N', 'B'],
            'D': ['D', 'B'],
            'Q': ['Q', 'Z'],
            'E': ['E', 'Z'],
            'J': ['I'],    
            'U': ['V']
        }    
                        
        new_search_nodes = []
        results = []
        while self.search_nodes[search_index]:
            current_node = self.search_nodes[search_index].pop() 
            #extra = []
            #if letter == 'N' or letter == 'D': extra = ['B']
            #if letter == 'Q' or letter == 'E': extra = ['Z']

            chars = ['X', 'O'] + special_chars.get(letter, [letter]) 
            #chars = [ letter ]
            #pset = list(filter(None, [ current_node.children[self._charToIndex(c)] for c in chars ]))
            pset = list(filter(None, [ current_node.children.get(c,None) for c in chars ]))
            #if not pset: 
            #    return ([], False)
            new_search_nodes += pset 
            results += [node.string for node in pset if node.string is not None ] 
            
        if not new_search_nodes: return ([], False)
        self.search_nodes[search_index] = new_search_nodes
        return (results, True) 
        #return pCrawl != None and pCrawl.isEndOfWord 



    def search(self, key): 
          
        # Search key in the trie 
        # Returns true if key presents  
        # in trie, else false 
        pCrawl = self.root 
        length = len(key) 
        for level in range(length): 
            index = self._charToIndex(key[level]) 
            if not pCrawl.children[index]: 
                return False
            pCrawl = pCrawl.children[index] 
  
        return pCrawl != None and pCrawl.isEndOfWord 
