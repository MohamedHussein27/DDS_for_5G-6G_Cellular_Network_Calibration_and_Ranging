class Transaction:
    def __init__(self, tid=0, data=0):
        self.tid = tid
        self.data = data

    def display(self):
        print("----- Transaction -----")
        print(f"TID           : {self.tid}")
        print(f"Data          : {self.data}")
        print("-----------------------")
