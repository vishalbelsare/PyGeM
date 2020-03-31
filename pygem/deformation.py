from abc import ABC, abstractmethod
 
class Deformation(ABC):
 
    @abstractmethod
    def __init__(self, value):
        pass
    
    @abstractmethod
    def __call__(self, src):
        pass
