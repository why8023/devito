class Boundary(Grid, BoundaryFunction):

    """
    An object that contains the data relevant for implementing the
    immersed boundary oject on the given domain.

    :param param0: Description.
    :param param1: (Optional) Description.

    Note: To add.

    """

    def __init__(self):
    
        # Step1: Check what kind of boundary function we have.
        # To start, this will only work for boundaries of the form
        # y_b = f(x).
        
        print(callable(BoundaryFunction))


    @property
    def method1(self):
        """My method."""
        return self._method1

 
