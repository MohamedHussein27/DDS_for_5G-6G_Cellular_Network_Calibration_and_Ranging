# transaction.py
class Transaction:
    def __init__(self, tid, data):
        self.tid = tid
        self.data = data
        self.processed_data = None