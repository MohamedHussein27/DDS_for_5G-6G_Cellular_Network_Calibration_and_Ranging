class Transaction:
    def __init__(self, tid, data):
        self.tid = tid
        self.data = data

    def __repr__(self):
        return f"Transaction(TID={self.tid}, DATA={self.data})"