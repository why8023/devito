from devito import Grid, Function

import numpy as np

import pandas as pd

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
        # To start, this will only work for boundaries of the form:
        # x_b = const : 1d.
        # y_b = f(x) : 2d.
        
        if not callable(BoundaryFunction):
            raise NotImplementedError
        
        # FIX ME: All this needs refactoring.
        self._primary_nodes(Grid, BoundaryFunction)
        
        # For 2D case:
        if np.asarray(Grid.shape).size == 1:
            # Add 1D case
            self._node_list = None
        elif np.asarray(Grid.shape).size == 2:
            # Now work out the full list
            self._node_list(Grid)
        else:
            raise NotImplementedError
            
        #print(self._secondary_nodes)
        

    def _primary_nodes(self, grid, BoundaryFunction):
        """Compute 'primary boundaru nodes."""
        
        # Should all this be done symbolically or not?

        if not np.ndim(grid.dimensions) <= 2:
            raise NotImplementedError
        
        dimensions = grid.dimensions

        # FIX ME: Must be a better way of doing this:        
        shape = np.asarray(grid.shape)
        extent = np.asarray(grid.extent)
        
        spacing = extent/(shape-1)
        
        if shape.size == 1:
            x_coords = np.linspace(0,extent[0],shape[0])
            # In this case the boundary is a single node
            boundary = BoundaryFunction()
            pn = np.floor(boundary/spacing[0]).astype(int)
            # Check if node is on or off the boundary:
            boundary_logic_list = None
            # FIX ME: Below code is horrible!!
            # First two cases shouldn't occur in 1D.
            if pn >= shape[0]:
                pass
            elif pn < 0:
                pass
            else:
                if abs(boundary-x_coords[pn]) < 2*np.pi*np.finfo(float).eps:
                    boundary_logic_list = 'on'
                else:
                    boundary_logic_list = 'off'
        elif shape.size == 2:
            # FIX ME: Support non zero origin case
            x_coords = np.linspace(0,extent[0],shape[0])
            y_coords = np.linspace(0,extent[1],shape[1])
            boundary = BoundaryFunction(x_coords)
            pn = np.floor(boundary/spacing[1]).astype(int)
            pn[pn < 0] = -1 # Replace this with some 'None' equiv
            pn[pn >= shape[1]] = -1
            # FIX ME: Disgusting code
            boundary_logic_list = [None]*shape[0]
            for j in range(0,shape[0]):
                if pn[j] >= shape[0]:
                    pass
                elif pn[j] < 0:
                    pass
                else:
                    if abs(boundary[j]-y_coords[pn[j]]) < 2*np.pi*np.finfo(float).eps:
                        boundary_logic_list[j] = 'on'
                    else:
                        boundary_logic_list[j] = 'off'
        else:
            raise NotImplementedError
        
        self._primary_nodes = pn
        self.boundary_logic = boundary_logic_list
        
        return self._primary_nodes
    
    
    def _node_list(self, grid):
        
        pn = self._primary_nodes
        
        dpn = np.zeros((pn.size,), dtype=int)
        
        for j in range(1,pn.size):
            if (pn[j] < 0) or (pn[j-1] < 0):
                dpn[j-1] = pn.size # Our 'int NaN'
            else:
                dpn[j-1] = pn[j]-pn[j-1]
        if dpn[-2] == pn.size:
            dpn[-1] = pn.size
        
        # dpn Cases:
        # 0: x stencil doesn't need modifying
        # 1: x stencil needs modifying on the left
        # 2+: x stencil +dpn-1 points below the primary points
        #     require left modification
        # -1+: As +ve case but with right modification instead.
        
        # Now build the full modification list + logic
        # Logic options:
        # 'y': Modify y stencil only
        # 'xl': Modify x stencil only on the left
        # 'xr': Modify x stencil only on the right
        # 'xm': Modify x stencil only on both left and right
        # 'yxl': Modify y stencil and x on the left
        # 'yxr': Modify y stencil and x on the right
        # 'yxm': Modify y stencil and x on both left and right
        
        count = 0
        
        # This isn't going to work:
        # We need to do the below via a recursive object
        # BETTER WAY: Create a grid (with redundancy) around the boundary
        # and just compute the eta's. A general stencil can then be worked
        # out with this info.
        coordinate_dict = ()
        coordinate_logic = ()
        for j in range(0,pn.size):
            # Add primary node + logic
            coordinate_dict = coordinate_dict + ([j,pn[j]],)
            if dpn[j] == pn.size:
                coordinate_logic = coordinate_logic + ('DoNothing',)
            elif dpn[j] == 0:
                coordinate_logic = coordinate_logic + ('y',)
            elif dpn[j] > 0:
                coordinate_logic = coordinate_logic + ('yxl',)
            elif dpn[j] < 0:
                coordinate_logic = coordinate_logic + ('yxr',)
            # Add other required nodes with this x-coordinate
            if coordinate_logic[count] == 'DoNothing':
                count+=1
            elif coordinate_logic[count] == 'y':
                coordinate_dict = coordinate_dict + ([j,pn[j]-1],)
                coordinate_logic = coordinate_logic + ('y',)
                count+=2
            elif coordinate_logic[count] == 'yxl':
                coordinate_dict = coordinate_dict + ([j,pn[j]-1],)
                coordinate_logic = coordinate_logic + ('y',)
                count+=2
            elif coordinate_logic[count] == 'yxr':
                coordinate_dict = coordinate_dict + ([j,pn[j]-1],)
                #coordinate_logic = coordinate_logic + ('y',)
                coordinate_logic = coordinate_logic + (['yxl','yxr'],)
                count+=2
            else:
                count+=1
        
        print(count)
            
        coordinate = pd.Series(coordinate_dict)
        logic = pd.Series(coordinate_logic)
        nodes = pd.DataFrame({'coordinate': coordinate,
                              'logic': logic})
        
        # So data structure should be
        # nodes = pd.DataFrame({'coordinate': coordinate,
        #                       'etax': etax, 'etay': etay})
        # where etay/etax can possibly be 'None' and etax can possibly
        # be multivalued (etay cannot be multivalued).
        
        print(nodes)
        
        self._secondary_nodes = nodes
        
        return self._secondary_nodes

 
