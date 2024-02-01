"""
The class file for global configuration.

Only a few settings can be changed here, and 
the parameters in this file are meant for a
temporary solution for experimental features.

Be cautious in adding more parameters here.
"""

class Config:

    __config = {
        "TIHamiltonian_tol" : None
    }

    @staticmethod
    def value(name):
        return Config.__config[name]
    
    @staticmethod
    def set(name, value):
        if name not in Config.__config:
            raise Exception("This configuration does not exist.")
        Config.__config[name] = value
