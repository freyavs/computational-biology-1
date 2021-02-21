import numpy as np

class TrieNode: 
      
    # Trie node class 
    def __init__(self): 
        self.children = [None]*26
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
  
    def _charToIndex(self,ch): 
          
        # private helper function 
        # Converts key current character into index 
        # use only 'a' through 'z' and lower case 
        #if ch == 'J':
            #return ord('I') - ord('A')
        #if ch == 'U': 
            #return ord('V') - ord('A')

        return ord(ch)-ord('A') 
  
  
    def insert(self,key): 
          
        # If not present, inserts key into trie 
        # If the key is prefix of trie node,  
        # just marks leaf node 
        pCrawl = self.root 
        length = len(key) 
        for level in range(length): 
            index = self._charToIndex(key[level]) 
  
            # if current character is not present 
            if not pCrawl.children[index]: 
                pCrawl.children[index] = self.getNode() 
            pCrawl = pCrawl.children[index] 
  
        # mark last node as leaf 
        pCrawl.string = key
        pCrawl.isEndOfWord = True

    def search_update(self, letter, search_index):
        pCrawl = self.search_nodes[search_index].pop() 
        index = self._charToIndex(letter) 
        if not pCrawl.children[index]: 
            return ([], False, False)
        pCrawl = pCrawl.children[index] 
        self.search_nodes[search_index].append(pCrawl)
        
        return ([ pCrawl.string ], pCrawl.isEndOfWord, True) 
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

# driver function 
def main(): 
  
    # Input keys (use only 'a' through 'z' and lower case) 
    keys = ["the","a","there","anaswe","any", 
            "by","their"] 
    output = ["Not present in trie", 
              "Present in trie"] 
  
    # Trie object 
    t = Trie() 
  
    # Construct trie 
    for key in keys: 
        t.insert(key) 
  
    # Search for different keys 
    print("{} ---- {}".format("the",output[t.search("the")])) 
    print("{} ---- {}".format("th",output[t.search("th")])) 
    print("{} ---- {}".format("these",output[t.search("these")])) 
    print("{} ---- {}".format("their",output[t.search("their")])) 
    print("{} ---- {}".format("thaw",output[t.search("thaw")])) 
  
if __name__ == '__main__': 
    main() 
  
# This code is contributed by Atul Kumar (www.facebook.com/atul.kr.007) 