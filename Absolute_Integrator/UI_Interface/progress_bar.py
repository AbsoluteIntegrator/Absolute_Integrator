from abc import ABCMeta, abstractmethod

class ProgressBarBase(metaclass=ABCMeta):
    @abstractmethod
    def set_title(self, title):
        """
        Sets informative text about what the progress bar is indicating the progress of
        """
        pass
    
    @abstractmethod
    def set_position(self, position):
        """
        Sets the current thing the progress bar is on
        """
        pass

    @abstractmethod
    def set_end(self, end):
        """
        Sets the total number of things for the progress bar to count to
        """
        pass
