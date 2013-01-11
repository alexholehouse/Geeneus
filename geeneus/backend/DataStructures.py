import collections

class CaseInsensitiveDict(collections.Mapping):
    def __init__(self, d):
        self._d = d
        self._s = dict((k.upper(), k) for k in d)
        self._lc = dict((k.upper(), d[k]) for k in d)
 
    def __contains__(self, k):
        return k.upper() in self._s
    
    def __len__(self):
        return len(self._s)
    
    def __iter__(self):
        return iter(self._s)
    
    def __getitem__(self, k):
        return self._d[self._s[k.upper()]]
    
    def actual_key_case(self, k):
        return self._s.get(k.upper())

    def __repr__(self): 
        return str(self._lc)
