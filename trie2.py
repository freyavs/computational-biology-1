import numpy as np

class TrieNode: 
      
    # Trie node class 
    def __init__(self): 
        self.children = [None]*26
  
        # isEndOfWord is True if node represent the end of the word 
        self.isEndOfWord = False
  
class Trie: 
      
    # Trie data structure class 
    def __init__(self): 
        self.root = self.getNode() 
  
    def getNode(self): 
      
        # Returns new trie node (initialized to NULLs) 
        return TrieNode() 
  
    def _charToIndex(self,ch): 
          
        # private helper function 
        # Converts key current character into index 
        # use only 'a' through 'z' and lower case 
        if ch == 'J':
            return ord('I') - ord('A')
        if ch == 'U': 
            return ord('V') - ord('A')

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
        pCrawl.isEndOfWord = True

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

    def search3(self, key): 
          
        # Search key in the trie 
        # Returns true if key presents  
        # in trie, else false 
        length = len(key) 
        zoekfront = [(self.root,0, "")]
        res = []

        while len(zoekfront) > 0:
            zoek = zoekfront.pop(0)
            pCrawl = zoek[0]
            current = zoek[2] # huidig pad van string

            for level in range(zoek[1], length): 
                chars = [key[level], 'O','X']
                pset = list(filter(None, [ pCrawl.children[self._charToIndex(c)] for c in chars ]))
                n = pCrawl.children[index] 
                if len(pset) == 0: 
                    break

                for i,p in enumerate(pset):
                    if i ==0: pCrawl = p
                    else: zoekfront.append(
    
            if pCrawl != None and pCrawl.isEndOfWord: res.append(current)

    def search2(self, key, startnode = None): 
        # Search key in the trie 
        # Returns true if key presents  
        # in trie, else false 
        if startnode: pCrawl = startnode
        else: pCrawl = self.root
        length = len(key) 
        for level in range(length): 
            index = self._charToIndex(key[level]) 
            #indices = [index, self._charToIndex('O'), self._charToIndex('X')] # o and x are wildcards and can be any letter in the protein we are searching
            indices = [index]
            indices_found = [ pCrawl.children[i] is not None for i in indices]
            indices_true = sum(indices_found)
            if indices_true == 0: return False
            indices_to_check = np.where(indices_found)[0]
            res = []

            if indices_true == 1: 
                pCrawl = pCrawl.children[indices[indices_to_check[0]]] 
                #return self.search(key[level+1:], pCrawl.children[index])

            else:
                for ind in indices_to_check:
                    res.append(self.search(key[level+1:], pCrawl.children[ind]))
                    return res
                
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