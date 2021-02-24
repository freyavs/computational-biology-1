import numpy as np

class Node: 
    def __init__(self): 
        self.children = dict()
        self.string = None 
  
class Trie: 
    def __init__(self): 
        self.special_chars = { 
            'N': ['B'],
            'D': ['B'],
            'Q': ['Z'],
            'E': ['Z'],
            'I': ['J'],    
            'V': ['U']
        }    
        self.root = Node() 
        self.search_nodes = [ [self.root], [self.root] ]

    def reset(self):
        self.search_nodes = [ [self.root], [self.root] ]
  
    def add(self,word): 
        current = self.root 
        for _, letter in enumerate(word): 
            if not current.children.get(letter, None): 
                current.children[letter] = Node() 
            current = current.children[letter] 
  
        current.string = word 

    def find_update(self, letter, search_index):
        new_search_nodes = []
        results = []
        search = self.search_nodes[search_index]
        chars = [letter, 'X', 'O']  + self.special_chars.get(letter, []) 

        while search:
            current_node = search.pop() 
            node_list = list(filter(None, [ current_node.children.get(c,None) for c in chars ]))
            new_search_nodes += node_list
            if node_list: results += list(filter(None, [node.string for node in node_list])) 
            
        if not new_search_nodes: return ([], False)
        self.search_nodes[search_index] = new_search_nodes

        return (results, True) 

    def find_update_simple(self, letter, search_index):
        current = self.search_nodes[search_index].pop() 
        if not current.children.get(letter, None): 
            return ([], False)
        current = current.children[letter] 
        self.search_nodes[search_index].append(current)
        
        return ([ current.string ], True) 

    def search(self, key): 
        pCrawl = self.root 
        length = len(key) 
        for level in range(length): 
            letter = key[level]
            if not pCrawl.children.get(letter, None): 
                return False
            pCrawl = pCrawl.children[letter] 
  
        return pCrawl != None and pCrawl.string 