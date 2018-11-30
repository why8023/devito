import weakref
import abc
import gc
from collections import namedtuple
from operator import mul
from functools import reduce
from ctypes import POINTER, Structure, byref

import numpy as np
import sympy
from cached_property import cached_property
from cgen import Struct, Value

from devito.data import DOMAIN, OWNED, HALO, LEFT, default_allocator
from devito.symbolics import Add
from devito.tools import (ArgProvider, EnrichedTuple, Pickable, ctypes_to_cstr,
                          dtype_to_cstr, dtype_to_ctype, memoized_meth)

__all__ = ['Symbol', 'Indexed']

# This cache stores a reference to each created data object
# so that we may re-create equivalent symbols during symbolic
# manipulation with the correct shapes, pointers, etc.
_SymbolCache = {}


class Basic(object):
    """
    Three relevant types inherit from this class: ::

        * AbstractSymbol: represents a scalar; may carry data; may be used
                          to build equations.
        * AbstractFunction: represents a discrete R^n -> R function; may
                            carry data; may be used to build equations.
        * AbstractObject: represents a generic object, for example a (pointer
                          to) data structure.

                                        Basic
                                          |
                    ------------------------------------------
                    |                     |                  |
             AbstractSymbol       AbstractFunction     AbstractObject

    Subtypes must implement a number of methods/properties to enable code
    generation via the Devito compiler. These methods/properties are easily
    recognizable as their name starts with _C_.

    Notes
    -----
    Part of the :class:`AbstractFunction` sub-hierarchy is implemented in
    :mod:`function.py`.
    """

    # Top hierarchy
    is_AbstractFunction = False
    is_AbstractSymbol = False
    is_AbstractObject = False

    # Symbolic objects created internally by Devito
    is_Symbol = False
    is_Array = False
    is_Object = False
    is_LocalObject = False

    # Created by the user
    is_Input = False
    # Scalar symbolic objects created by the user
    is_Dimension = False
    is_Constant = False
    # Tensor symbolic objects created by the user
    is_TensorFunction = False
    is_Function = False
    is_TimeFunction = False
    is_SparseTimeFunction = False
    is_SparseFunction = False
    is_PrecomputedSparseFunction = False
    is_PrecomputedSparseTimeFunction = False

    # Basic symbolic object properties
    is_Scalar = False
    is_Tensor = False

    @abc.abstractmethod
    def __init__(self, *args, **kwargs):
        return

    @abc.abstractproperty
    def _C_name(self):
        """
        The C-level name of the object.

        Returns
        -------
        str
        """
        return

    @abc.abstractproperty
    def _C_typename(self):
        """
        The C-level type of the object.

        Returns
        -------
        str
        """
        return

    @abc.abstractproperty
    def _C_typedata(self):
        """
        The C-level type of the data values.

        Returns
        -------
        str
        """
        return

    @abc.abstractproperty
    def _C_ctype(self):
        """
        The C-level type of the object, as a ctypes object, suitable for type
        checking when calling functions via ctypes.

        Returns
        -------
        ctypes type
        """
        return

    @property
    def _C_typedecl(self):
        """
        The C-level struct declaration representing the object.

        Returns
        -------
        :class:`cgen.Struct` or None
            None if the object C type can be expressed with a basic C type,
            such as float or int.
        """
        return


class Cached(object):
    """
    Base class for symbolic objects that cache on the class type.

    In order to maintain meta information across the numerous
    re-instantiation SymPy performs during symbolic manipulation, we inject
    the symbol name as the class name and cache all created objects on that
    name. This entails that a symbolic object inheriting from :class:`Cached`
    should implement `__init__` in the following way:

        .. code-block::
            def __init__(self, \*args, \*\*kwargs):
                if not self._cached():
                    ... # Initialise object properties from kwargs
    """

    @classmethod
    def _cached(cls):
        """Test if current class is already in the symbol cache."""
        return cls in _SymbolCache

    @classmethod
    def _cache_put(cls, obj):
        """Store given object instance in symbol cache.

        :param obj: Object to be cached.
        """
        _SymbolCache[cls] = weakref.ref(obj)

    @classmethod
    def _symbol_type(cls, name):
        """Create new type instance from cls and inject symbol name"""
        return type(name, (cls, ), dict(cls.__dict__))

    def _cached_init(self):
        """Initialise symbolic object with a cached object state"""
        original = _SymbolCache[self.__class__]
        self.__dict__ = original().__dict__

    def __hash__(self):
        """The hash value of an object that caches on its type is the
        hash value of the type itself."""
        return hash(type(self))


class AbstractSymbol(sympy.Symbol, Basic, Pickable):
    """
    Base class for dimension-free symbols, only cached by SymPy.

    The sub-hierarchy is structured as follows

                             AbstractSymbol
                                   |
                 -------------------------------------
                 |                                   |
        AbstractCachedSymbol                      Dimension
                 |
        --------------------
        |                  |
     Symbol            Constant
        |
     Scalar

    There are three relevant :class:`AbstractSymbol` sub-types: ::

        * Symbol: A generic scalar symbol that can be used to build an equation.
                  It does not carry data. Typically, :class:`Symbol`s are
                  created internally by Devito (e.g., for temporary variables)
        * Constant: A generic scalar symbol that can be used to build an equation.
                    It carries data (a scalar value).
        * Dimension: A problem dimension, used to create an iteration space. It
                     may be used to build equations; typically, it is used as
                     an index for a :class:`Indexed`.
    """

    is_AbstractSymbol = True

    @property
    def indices(self):
        return ()

    @property
    def shape(self):
        return ()

    @property
    def ndim(self):
        return 0

    @property
    def symbolic_shape(self):
        return ()

    @property
    def base(self):
        return self

    @property
    def function(self):
        return self

    def indexify(self):
        return self

    def _subs(self, old, new, **hints):
        """This stub allows sympy.Basic.subs to operate on an expression
        involving devito Scalars.  Ordinarily the comparisons between
        devito subclasses of sympy types are quite strict."""
        try:
            if old.name == self.name:
                return new
        except AttributeError:
            pass

        return self


class AbstractCachedSymbol(AbstractSymbol, Cached):
    """
    Base class for dimension-free symbols, cached by both Devito and SymPy.

    For more information, refer to the documentation of :class:`AbstractSymbol`.
    """

    def __new__(cls, *args, **kwargs):
        options = kwargs.get('options', {})
        if cls in _SymbolCache:
            newobj = sympy.Symbol.__new__(cls, *args, **options)
            newobj._cached_init()
        else:
            name = kwargs.get('name')

            # Create the new Function object and invoke __init__
            newcls = cls._symbol_type(name)
            newobj = sympy.Symbol.__new__(newcls, name, *args, **options)

            # Initialization
            newobj._dtype = cls.__dtype_setup__(**kwargs)
            newobj.__init__(*args, **kwargs)

            # Store new instance in symbol cache
            newcls._cache_put(newobj)
        return newobj

    __hash__ = Cached.__hash__

    @classmethod
    def __dtype_setup__(cls, **kwargs):
        """Extract the object data type from ``kwargs``."""
        return None

    @property
    def dtype(self):
        """The data type of the object."""
        return self._dtype

    @property
    def _C_name(self):
        return self.name

    @property
    def _C_typename(self):
        return 'const %s' % dtype_to_cstr(self.dtype)

    @property
    def _C_typedata(self):
        return dtype_to_cstr(self.dtype)

    @property
    def _C_ctype(self):
        return dtype_to_ctype(self.dtype)

    # Pickling support
    _pickle_kwargs = ['name', 'dtype']
    __reduce_ex__ = Pickable.__reduce_ex__

    @property
    def _pickle_reconstruct(self):
        return self.__class__.__base__


class Symbol(AbstractCachedSymbol):

    """A :class:`sympy.Symbol` capable of mimicking an :class:`sympy.Indexed`"""

    is_Symbol = True

    @classmethod
    def __dtype_setup__(cls, **kwargs):
        return kwargs.get('dtype', np.int32)


class Scalar(Symbol):
    """Symbolic object representing a scalar.

    :param name: Name of the symbol.
    :param dtype: (Optional) data type of the object. Defaults to float32.
    """

    is_Scalar = True

    @classmethod
    def __dtype_setup__(cls, **kwargs):
        return kwargs.get('dtype', np.float32)

    @property
    def _mem_stack(self):
        """Return True if the associated data should be allocated on the stack
        in a C module, False otherwise."""
        return True

    def update(self, dtype=None, **kwargs):
        self.dtype = dtype or self.dtype


class AbstractFunction(sympy.Function, Basic, Pickable):
    """
    Base class for tensor symbols, only cached by SymPy. It inherits from and
    mimick the behaviour of a :class:`sympy.Function`.

    The sub-hierarchy is structured as follows

                         AbstractFunction
                                |
                      AbstractCachedFunction
                                |
                 ---------------------------------
                 |                               |
           TensorFunction                      Array
                 |
         ----------------------------------------
         |                                      |
         |                           AbstractSparseFunction
         |                                      |
         |               -----------------------------------------------------
         |               |                      |                            |
      Function     SparseFunction   AbstractSparseTimeFunction  PrecomputedSparseFunction
         |               |                      |                            |
         |               |   ------------------------------------     --------
         |               |   |                                  |     |
    TimeFunction  SparseTimeFunction                 PrecomputedSparseTimeFunction

    There are five relevant :class:`AbstractFunction` sub-types: ::

        * Array: A function that does not carry data. Usually created by the DSE.
        * Function: A space-varying discrete function, which carries user data.
        * TimeFunction: A time- and space-varying discrete function, which carries
                        user data.
        * SparseFunction: A space-varying discrete function representing "sparse"
                          points, i.e. points that are not aligned with the
                          computational grid.
        * SparseTimeFunction: A time- and space-varying function representing "sparse"
                          points, i.e. points that are not aligned with the
                          computational grid.
        * PrecomputedSparseFunction: A SparseFunction that uses a custom interpolation
                                     scheme, instead of the included linear interpolators.
        * PrecomputedSparseTimeFunction: A SparseTimeFunction that uses a custom
                                         interpolation scheme, instead of the included
                                         linear interpolators.
    """

    is_AbstractFunction = True


class AbstractCachedFunction(AbstractFunction, Cached):
    """
    Base class for tensor symbols, cached by both Devito and Sympy.

    For more information, refer to the documentation of :class:`AbstractFunction`.
    """

    def __new__(cls, *args, **kwargs):
        options = kwargs.get('options', {})
        if cls in _SymbolCache:
            newobj = sympy.Function.__new__(cls, *args, **options)
            newobj._cached_init()
        else:
            name = kwargs.get('name')
            indices = cls.__indices_setup__(**kwargs)

            # Create the new Function object and invoke __init__
            newcls = cls._symbol_type(name)
            newobj = sympy.Function.__new__(newcls, *indices, **options)

            # Initialization. The following attributes must be available
            # when executing __init__
            newobj._name = name
            newobj._indices = indices
            newobj._shape = cls.__shape_setup__(**kwargs)
            newobj._dtype = cls.__dtype_setup__(**kwargs)
            newobj.__init__(*args, **kwargs)

            # All objects cached on the AbstractFunction /newobj/ keep a reference
            # to /newobj/ through the /function/ field. Thus, all indexified
            # object will point to /newobj/, the "actual Function".
            newobj.function = newobj

            # Store new instance in symbol cache
            newcls._cache_put(newobj)
        return newobj

    def __init__(self, *args, **kwargs):
        if not self._cached():
            # Setup halo and padding regions
            self._is_halo_dirty = False
            self._in_flight = []
            self._halo = self.__halo_setup__(**kwargs)
            self._padding = self.__padding_setup__(**kwargs)

    __hash__ = Cached.__hash__

    @classmethod
    def __indices_setup__(cls, **kwargs):
        """Extract the object indices from ``kwargs``."""
        return ()

    @classmethod
    def __shape_setup__(cls, **kwargs):
        """Extract the object shape from ``kwargs``."""
        return ()

    @classmethod
    def __dtype_setup__(cls, **kwargs):
        """Extract the object data type from ``kwargs``."""
        return None

    def __halo_setup__(self, **kwargs):
        return kwargs.get('halo', tuple((0, 0) for i in range(self.ndim)))

    def __padding_setup__(self, **kwargs):
        return kwargs.get('padding', tuple((0, 0) for i in range(self.ndim)))

    @property
    def name(self):
        """Return the name of the object."""
        return self._name

    @property
    def indices(self):
        """Return the indices (aka dimensions) of the object."""
        return self._indices

    @property
    def dimensions(self):
        """Tuple of :class:`Dimension`s representing the object indices."""
        return self.indices

    @property
    def shape(self):
        """Return the shape of the object."""
        return self._shape

    @property
    def dtype(self):
        """Return the data type of the object."""
        return self._dtype

    @property
    def ndim(self):
        """Return the rank of the function."""
        return len(self.indices)

    @property
    def symbolic_shape(self):
        """
        The symbolic shape of the object. This includes the domain, halo, and
        padding regions. While halo and padding are known quantities (integers),
        the domain size is represented by a symbol.
        """
        halo_sizes = [Add(*i) for i in self._extent_halo]
        pad_sizes = [Add(*i) for i in self._extent_padding]
        domain_sizes = [i.symbolic_size for i in self.indices]
        ret = tuple(Add(i, j, k) for i, j, k in zip(domain_sizes, halo_sizes, pad_sizes))
        return EnrichedTuple(*ret, getters=self.dimensions)

    @property
    def indexed(self):
        """Extract a :class:`IndexedData` object from the current object."""
        return IndexedData(self.name, shape=self.shape, function=self.function)

    @property
    def _mem_external(self):
        """Return True if the associated data was/is/will be allocated directly
        from Python (e.g., via NumPy arrays), False otherwise."""
        return False

    @property
    def _mem_stack(self):
        """Return True if the associated data was/is/will be allocated on the stack
        in a C module, False otherwise."""
        return False

    @property
    def _mem_heap(self):
        """Return True if the associated data was/is/will be allocated on the heap
        in a C module, False otherwise."""
        return False

    @property
    def size(self):
        """Return the number of elements this function is expected to store in memory.
           Note that this would need to be combined with self.dtype to give the actual
           size in bytes
        """
        return reduce(mul, self.shape)

    @property
    def halo(self):
        return self._halo

    @property
    def padding(self):
        return self._padding

    @property
    def _C_name(self):
        return "%s_vec" % self.name

    @property
    def _C_typedata(self):
        return dtype_to_cstr(self.dtype)

    @cached_property
    def _offset_domain(self):
        """
        The number of points between the first (last) allocated element
        and the first (last) domain element, for each dimension.
        """
        left = tuple(np.add(self._extent_halo.left, self._extent_padding.left))
        right = tuple(np.add(self._extent_halo.right, self._extent_padding.right))

        Offset = namedtuple('Offset', 'left right')
        offsets = tuple(Offset(i, j) for i, j in np.add(self._halo, self._padding))

        return EnrichedTuple(*offsets, getters=self.dimensions, left=left, right=right)

    @cached_property
    def _offset_halo(self):
        """
        The number of points between the first (last) allocated element
        and the first (last) halo element, for each dimension.
        """
        left = self._extent_padding.left
        right = self._extent_padding.right

        Offset = namedtuple('Offset', 'left right')
        offsets = tuple(Offset(i, j) for i, j in self._padding)

        return EnrichedTuple(*offsets, getters=self.dimensions, left=left, right=right)

    @cached_property
    def _extent_domain(self):
        """
        The number of points in the domain region.
        """
        return EnrichedTuple(*self.shape, getters=self.dimensions)

    @cached_property
    def _extent_halo(self):
        """
        The number of points in the halo region.
        """
        left = tuple(zip(*self._halo))[0]
        right = tuple(zip(*self._halo))[1]

        Extent = namedtuple('Extent', 'left right')
        extents = tuple(Extent(i, j) for i, j in self._halo)

        return EnrichedTuple(*extents, getters=self.dimensions, left=left, right=right)

    @cached_property
    def _extent_padding(self):
        """
        The number of points in the padding region.
        """
        left = tuple(zip(*self._padding))[0]
        right = tuple(zip(*self._padding))[1]

        Extent = namedtuple('Extent', 'left right')
        extents = tuple(Extent(i, j) for i, j in self._padding)

        return EnrichedTuple(*extents, getters=self.dimensions, left=left, right=right)

    @memoized_meth
    def _region_meta(self, region, dim, side=None):
        """
        Offset, relative to the origin, and extent of a given data region.

        Parameters
        ----------
        region : DataRegion
            The data region of interest (e.g., OWNED, HALO).
        dim : Dimension
            The dimension of interest.
        side : DataSide, optional
            The side of interest (LEFT, RIGHT). Irrelevant for certain regions.
        """
        if region is DOMAIN:
            offset = self._offset_domain[dim].left
            extent = self._extent_domain[dim]
        elif region is OWNED:
            if side is LEFT:
                offset = self._offset_domain[dim].left
                extent = self._extent_halo[dim].right
            else:
                offset = self._offset_halo[dim].left + self._extent_domain[dim]
                extent = self._extent_halo[dim].left
        elif region is HALO:
            if side is LEFT:
                offset = self._offset_halo[dim].left
                extent = self._extent_halo[dim].left
            else:
                offset = self._offset_domain[dim].left + self._extent_domain[dim]
                extent = self._extent_halo[dim].right
        else:
            raise ValueError("Unknown region `%s`" % str(region))

        RegionMeta = namedtuple('RegionMeta', 'offset extent')
        return RegionMeta(offset, extent)

    @property
    def _data_alignment(self):
        """The address of any allocated memory is a multiple of the alignment."""
        return default_allocator().guaranteed_alignment

    def indexify(self, indices=None):
        """Create a :class:`sympy.Indexed` object from the current object."""
        if indices is not None:
            return Indexed(self.indexed, *indices)

        # Get spacing symbols for replacement
        spacings = [i.spacing for i in self.indices]

        # Only keep the ones used as indices.
        spacings = [s for i, s in enumerate(spacings)
                    if s.free_symbols.intersection(self.args[i].free_symbols)]

        # Substitution for each index
        subs = dict([(s, 1) for s in spacings])

        # Indices after substitutions
        indices = [a.subs(subs) for a in self.args]

        return Indexed(self.indexed, *indices)

    def __getitem__(self, index):
        """Shortcut for ``self.indexed[index]``."""
        return self.indexed[index]

    # Pickling support
    _pickle_kwargs = ['name', 'dtype', 'halo', 'padding']
    __reduce_ex__ = Pickable.__reduce_ex__

    @property
    def _pickle_reconstruct(self):
        return self.__class__.__base__


class Array(AbstractCachedFunction):
    """A symbolic function object, created and managed directly by Devito..

    :param name: Name of the object.
    :param dtype: Data type of the object.
    :param dimensions: The symbolic dimensions of the object.
    :param halo: The halo region of the object, expressed as an iterable
                 ``[(dim1_left_halo, dim1_right_halo), (dim2_left_halo, ...)]``
    :param padding: The padding region of the object, expressed as an iterable
                    ``[(dim1_left_pad, dim1_right_pad), (dim2_left_pad, ...)]``
    :param scope: (Optional) Control memory allocation. Allowed values are
                  ['heap', 'stack', 'external']. Defaults to 'heap'.
                  'external' implies that no storage needs to be allocated
                  for the Array.
    """

    is_Array = True
    is_Tensor = True

    def __new__(cls, *args, **kwargs):
        kwargs.update({'options': {'evaluate': False}})
        return AbstractCachedFunction.__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        if not self._cached():
            super(Array, self).__init__(*args, **kwargs)

            self._scope = kwargs.get('scope', 'heap')
            assert self._scope in ['heap', 'stack', 'external']

    @classmethod
    def __indices_setup__(cls, **kwargs):
        return tuple(kwargs['dimensions'])

    @classmethod
    def __dtype_setup__(cls, **kwargs):
        return kwargs.get('dtype', np.float32)

    @property
    def shape(self):
        return self.symbolic_shape

    @property
    def _mem_external(self):
        return self._scope == 'external'

    @property
    def _mem_stack(self):
        return self._scope == 'stack'

    @property
    def _mem_heap(self):
        return self._scope == 'heap'

    @property
    def _C_typename(self):
        return ctypes_to_cstr(POINTER(dtype_to_ctype(self.dtype)))

    def update(self, **kwargs):
        self._shape = kwargs.get('shape', self.shape)
        self._indices = kwargs.get('dimensions', self.indices)
        self._dtype = kwargs.get('dtype', self.dtype)
        self._halo = kwargs.get('halo', self._halo)
        self._padding = kwargs.get('padding', self._padding)
        self._scope = kwargs.get('scope', self._scope)
        assert self._scope in ['heap', 'stack', 'external']

    _pickle_kwargs = ['name', 'halo', 'padding', 'dimensions']


# Objects belonging to the Devito API not involving data, such as data structures
# that need to be passed to external libraries


class AbstractObject(Basic, sympy.Basic, Pickable):

    """
    Represent a generic pointer object.
    """

    is_AbstractObject = True

    def __new__(cls, *args, **kwargs):
        obj = sympy.Basic.__new__(cls)
        obj.__init__(*args, **kwargs)
        return obj

    def __init__(self, name, dtype):
        self.name = name
        self.dtype = dtype

    def __repr__(self):
        return self.name

    __str__ = __repr__

    def _hashable_content(self):
        return (self.name, self.dtype)

    @property
    def free_symbols(self):
        return {self}

    @property
    def _C_name(self):
        return self.name

    @property
    def _C_typename(self):
        return ctypes_to_cstr(self.dtype)

    @property
    def _C_ctype(self):
        return self.dtype

    # Pickling support
    _pickle_args = ['name', 'dtype']
    __reduce_ex__ = Pickable.__reduce_ex__


class Object(AbstractObject, ArgProvider):

    """
    Represent a generic pointer object, provided by the outside world.
    """

    is_Object = True

    def __init__(self, name, dtype, value=None):
        super(Object, self).__init__(name, dtype)
        self.value = value

    def _arg_defaults(self):
        if callable(self.value):
            return {self.name: self.value()}
        else:
            return {self.name: self.value}

    def _arg_values(self, **kwargs):
        if self.name in kwargs:
            return {self.name: kwargs.pop(self.name)}
        else:
            return {}


class CompositeObject(Object):

    """
    Represent a pointer object to a composite type (i.e., a C struct),
    provided by the outside world.
    """

    _dtype_cache = {}

    @classmethod
    def _generate_unique_dtype(cls, pname, pfields):
        dtype = POINTER(type(pname, (Structure,), {'_fields_': pfields}))
        key = (pname, tuple(pfields))
        return cls._dtype_cache.setdefault(key, dtype)

    def __init__(self, name, pname, pfields, value=None):
        dtype = CompositeObject._generate_unique_dtype(pname, pfields)
        value = value or byref(dtype._type_())
        super(CompositeObject, self).__init__(name, dtype, value)

    @property
    def pfields(self):
        return tuple(self.dtype._type_._fields_)

    @property
    def pname(self):
        return self.dtype._type_.__name__

    @property
    def fields(self):
        return [i for i, _ in self.pfields]

    def _hashable_content(self):
        return (self.name, self.pfields)

    @cached_property
    def _C_typedecl(self):
        return Struct(self.pname, [Value(ctypes_to_cstr(j), i) for i, j in self.pfields])

    # Pickling support
    _pickle_args = ['name', 'pname', 'pfields']
    _pickle_kwargs = []


class LocalObject(AbstractObject):

    """
    Represent a generic, yet locally defined, pointer object.
    """

    is_LocalObject = True


# Extended SymPy hierarchy follows, for essentially two reasons:
# - To keep track of `function`
# - To override SymPy caching behaviour


class IndexedData(sympy.IndexedBase, Pickable):
    """Wrapper class that inserts a pointer to the symbolic data object"""

    def __new__(cls, label, shape=None, function=None):
        obj = sympy.IndexedBase.__new__(cls, label, shape)
        obj.function = function
        return obj

    def func(self, *args):
        obj = super(IndexedData, self).func(*args)
        obj.function = self.function
        return obj

    def __getitem__(self, indices, **kwargs):
        """
        Return :class:`Indexed`, rather than :class:`sympy.Indexed`.
        """
        indexed = super(IndexedData, self).__getitem__(indices, **kwargs)
        return Indexed(*indexed.args)

    # Pickling support
    _pickle_kwargs = ['label', 'shape', 'function']
    __reduce_ex__ = Pickable.__reduce_ex__


class Indexed(sympy.Indexed):

    # The two type flags have changed in upstream sympy as of version 1.1,
    # but the below interpretation is used throughout the DSE and DLE to
    # identify Indexed objects. With the sympy-1.1 changes a new flag
    # obj.is_Indexed was introduced which should be preferred, but the
    # required changes are cumbersome and many...
    is_Symbol = False
    is_Atom = False

    def _hashable_content(self):
        return super(Indexed, self)._hashable_content() + (self.base.function,)

    @property
    def function(self):
        return self.base.function

    @property
    def dtype(self):
        return self.function.dtype

    @property
    def name(self):
        return self.function.name

# Utilities


class CacheManager(object):

    """
    Drop unreferenced objects from the SymPy and Devito caches. The associated
    data is lost (and thus memory is freed).
    """

    @classmethod
    def clear(cls):
        sympy.cache.clear_cache()
        gc.collect()
        for key, val in list(_SymbolCache.items()):
            if val() is None:
                del _SymbolCache[key]
