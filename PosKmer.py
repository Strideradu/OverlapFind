"""
class for kmer with position
"""
class PosKmer(object):
    def __init__(self, str, pos):
        self.word  = str
        self.pos = pos