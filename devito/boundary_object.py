from devito import Grid, Function

import numpy as np

class Boundary(object):

    """
    An object that contains the data relevant for implementing the
    immersed boundary oject on the given domain.

    :param param0: Description.
    :param param1: (Optional) Description.

    Note: To add.

    """

    def __init__(self, Grid, BoundaryFunction):
    
        # Step1: Check what kind of boundary function we have.
        # To start, this will only work for boundaries of the form
        # y_b = f(x).
        
        if not callable(BoundaryFunction):
            raise NotImplementedError
        
        self._primary_nodes(Grid, BoundaryFunction)
        

    def _primary_nodes(self, grid, BoundaryFunction):
        """Compute 'primary boundaru nodes."""
        
        # Should all this be done symbolically or not?

        if not np.ndim(grid.dimensions) <= 2:
            raise NotImplementedError
        
        dimensions = grid.dimensions

        # FIX ME: Must be a better way of doing this:        
        shape = grid.shape
        extent = grid.extent
        
        
        
        self._primary_nodes = None
        
        return self._primary_nodes

 
