class BaseProvider:
    def __init__(self):
        self.prog = None
        self.task = None
        self.qs_names = None

    @staticmethod
    def to_bin(res, n):
        b = bin(res)[2:]
        return "0" * (n - len(b)) + b

    def print_sites(self):
        if self.prog is None:
            raise Exception("No compiled job in record.")
        print("Order of sites:", self.qs_names)
