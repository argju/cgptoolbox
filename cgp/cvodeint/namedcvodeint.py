"""Cvode wrapper with named state variables and parameters."""

from contextlib import contextmanager
from collections import namedtuple

import numpy as np

from cvodeint.core import Cvodeint
from utils.dotdict import Dotdict

class Namedcvodeint(Cvodeint):
    """
    Cvode wrapper with named state variables and parameters.
    
    Constructor arguments are as for :class:`~cvodeint.core.Cvodeint`, 
    except that *p* is a recarray. Further ``*args, **kwargs``
    are passed to :class:`~cvodeint.core.Cvodeint`.
        
    With no arguments, this returns the van der Pol model as an example.
        
    >>> Namedcvodeint()
    Namedcvodeint(f_ode=vanderpol, t=array([0, 1]), y=[-2.0, 0.0])
        
    .. todo:: 
        
        Should be::
    
            Namedcvodeint(f_ode=vanderpol, t=array([0, 1]), 
            y=rec.array([(-2.0, 0.0)], dtype=[('x', '<f8'), ('y', '<f8')]), 
            p=rec.array([(1.0,)], dtype=[('epsilon', '<f8')]))

    CVODE keeps state variables in a structure of class NVector, which does not 
    offer named fields. The ODE right-hand side must therefore refer to state 
    variables by index, or convert from/to a record array.
    
    .. todo:: The cleanest design would be to share memory between a record array 
       and the NVector, but I haven't got that to work.
    
    This example is similar to the one for :class:`~cvodeint.core.Cvodeint`. 
    Call the :class:`Namedcvodeint` constructor with arguments 
    ODE function, time, initial state and parameter array.
    
    >>> ode, t, y, p = Namedcvodeint.example()
    >>> n = Namedcvodeint(ode, t, y, p)
    
    The resulting object has a :class:`Recarraylink` to the state NVector.
    
    >>> n.yr.x
    array([-2.])
    
    Parameters can be read or written via the record array.
    
    >>> n.pr.epsilon
    array([ 1.])
    """
    
    @staticmethod
    def example():
        """
        Example for :class:`Namedcvodeint`: van der Pol model.
        
        >>> Namedcvodeint.example()
        Example(ode=<function vanderpol at 0x...>, t=[0, 1], 
        y=rec.array([(-2.0, 0.0)], dtype=[('x', '<f8'), ('y', '<f8')]), 
        p=rec.array([(1.0,)], dtype=[('epsilon', '<f8')]))
        """
        t = [0, 1]
        y = np.rec.fromrecords([(-2.0, 0.0)], names="x y".split())
        p = np.rec.fromrecords([(1.0,)], names="epsilon")
        
        def vanderpol(t, y, ydot, f_data):
            ydot[0] = y[1]
            ydot[1] = p.epsilon * (1 - y[0] * y[0]) * y[1] - y[0]
        
        Nt = namedtuple("Example", "ode t y p")
        return Nt(vanderpol, t, y, p)
    
    def __init__(self, f_ode=None, t=None, y=None, p=None, *args, **kwargs):
        if f_ode is None:
            self.__init__(*self.example())
            return
        
        super(Namedcvodeint, self).__init__(f_ode, t, y.view(float), 
            *args, **kwargs)
        self.yr = Recarraylink(self.y, y.dtype)
        self.pr = p
        # objects that should not be reassigned, but whose value may change
        self.reassignwarning = """
        To change initial values or parameters for a model, 
        use model.y[:] = ... or model.pr[:] = ..., 
        not just model.y = ... or model.pr = ...
        The latter will cause the name model.y to point to a new object, 
        breaking the link to the CVODE object.
        With Numpy arrays, using x[:] = ... is guaranteed to modify 
        contents only.
        """
        self.originals = dict(pr=self.pr, y=self.y, yr=self.yr)
        self.dtype = Dotdict(y=y.dtype, p=p.dtype)
    
    def integrate(self, **kwargs):
        """
        Return Cvodeint.integrate() of CellML model; convert state to recarray

        :parameters: See :meth:`cvodeint.core.Cvodeint.integrate`
        :return tuple:
            * **t** : time vector
            * **Yr** : state recarray. Yr[i] is state at time t[i]
              Yr.V is state variable V; Yr.V[i] is V at time t[i]
            * **flag** : last flag returned by CVode
        
        To convert the recarray to a normal array, use ``Y = Yr.view(float)``.

        >>> vdp = Namedcvodeint() # default van der pol model
        >>> t, Yr, flag = vdp.integrate(t=np.linspace(0, 20, 100))
        >>> Y = Yr.view(float)
        >>> Yr[0]
        rec.array([(-2.0, 0.0)], dtype=[('x', '<f8'), ('y', '<f8')])        
        >>> Y[0]
        array([-2.,  0.])
        >>> Yr.x
        array([[-2.        ], [-1.96634283], ... [-1.96940322], [-2.00814991]])
        """
        if not all(self.__dict__[k] is v for k, v in self.originals.items()):
            raise AssertionError(self.reassignwarning)
        t, Y, flag = super(Namedcvodeint, self).integrate(**kwargs)
        Yr = Y.view(self.dtype.y, np.recarray)
        return t, Yr, flag
    
    @contextmanager
    def autorestore(self, _p=None, _y=None, **kwargs):
        """
        Context manager to restore state and parameters after use.
        
        :param array_like _p: Temporary parameter vector
        :param array_like _y: Temporary initial state
        :param dict ``**kwargs``: dict of (name, value) for parameters or initial state
        
        In the following example, the assignment to epsilon and the call 
        to :meth:`~cvodeint.core.Cvodeint.integrate`
        change the parameters and state of the model object.
        This change is undone by the :meth:`autorestore` context manager.
        
        >>> vdp = Namedcvodeint()
        >>> before = np.array([vdp.yr.x, vdp.pr.epsilon])
        >>> with vdp.autorestore():
        ...     vdp.pr.epsilon = 2
        ...     t, y, flag = vdp.integrate()
        >>> after = np.array([vdp.yr.x, vdp.pr.epsilon])
        >>> all(before == after)
        True
        
        .. note:: Rootfinding settings cannot be restored, and so are cleared
           on exiting the .autorestore() context manager
        
        Optionally, you can specify initial parameters and state as _p and _y. 
        Any key=value pairs are used to update parameters or state if the key 
        exists in exactly one of them, otherwise an error is raised.
        
        >>> pr = np.copy(vdp.pr).view(np.recarray)
        >>> y0 = np.copy(vdp.y).view(vdp.dtype.y, np.recarray)
        >>> pr.epsilon = 123
        >>> y0.x = 456
        >>> with vdp.autorestore(pr, y0, y=789):
        ...     vdp.pr.epsilon, vdp.yr.x, vdp.yr.y
        (array([ 123.]), array([ 456.]), array([ 789.]))
        >>> vdp.pr.epsilon, vdp.yr.x, vdp.yr.y
        (array([ 1.]), array([-2.]), array([ 0.]))
        
        .. todo:: Refactor pr and y into properties/OrderedDict objects with 
           update() methods.
        
        Settings are restored even in case of exception.
        
        >>> with vdp.autorestore():
        ...     vdp.pr.epsilon = 50
        ...     raise Exception
        Traceback (most recent call last):
        Exception
        >>> bool(vdp.pr.epsilon == before[1])
        True
        """
        oldy, oldpar = np.copy(self.y), np.copy(self.pr)
        if _p is not None:
            self.pr[:] = _p
        if _y is not None:
            try:
                self.y[:] = _y
            except ValueError: # can only convert an array of size 1 to a Python scalar
                self.y[:] = _y.squeeze()
            except TypeError: # float expected instead of numpy.void instance
                self.y[:] = _y.item()
        for k, v in kwargs.items():
            if k in self.dtype.p.names and k not in self.dtype.y.names:
                self.pr[k] = v
                continue
            if k in self.dtype.y.names and k not in self.dtype.p.names:
                # Recarraylink does not support item assignment
                setattr(self.yr, k, v)
                continue
            if k not in [self.dtype.y.names + self.dtype.p.names]:
                raise TypeError("Key %s not in parameter or rate vectors" % k)
            raise TypeError(
                "Key %s occurs in both parameter and state vectors" % k)
        
        try:
            yield
        finally:
            self.RootInit(0) # Disable any rootfinding
            self.y[:], self.pr[:] = oldy, oldpar
    
    # TODO: Enable dynamic regulation of selected parameters, eg stim_amplitude.
    # Create a new Namedcvodeint with a callback function as an attribute and 
    # ODE right-hand side as follows:
    #   Store original stim_amplitude
    #   Compute new stim_amplitude from t, y, p using callback function
    #   Compute ydot based on new stim_amplitude
    #   Restore original stim_amplitude
    #   Return ydot
    # Now the callback function can be substituted and integration should still 
    # work. The right-hand side could probably be included in the Cython code.
    # For now, just combine stim_amplitude etc (for pacing) and clamp/dynclamp.
    
    @contextmanager
    def clamp(self, **kwargs):
        """
        Derived model with state value(s) clamped to constant values.
        
        Names and values of clamped variables are given as keyword arguments.
        Clamped state variables cannot be modified along the way.
        
        Use as a context manager.
        
        >>> vdp = Namedcvodeint()
        >>> with vdp.clamp(x=0.5) as clamped:
        ...     t, y, flag = clamped.integrate(t=[0, 10])
        
        The clamped model has its own instance of the CVODE integrator and 
        state NVector.
        
        >>> clamped.cvode_mem is vdp.cvode_mem
        False
        
        However, changes in state are copied to the original on exiting the 
        'with' block.
        
        >>> vdp.yr.x
        array([ 0.5])
        
        The clamped model has the same parameter array as the original.
        You will usually not change parameters inside the 'with' block, 
        but if you do, it will affect both models. Be careful.
        
        >>> with vdp.clamp(x=0.5) as clamped:
        ...     clamped.pr.epsilon = 2
        >>> vdp.pr.epsilon
        array([ 2.])
        """
        # Indices to state variables whose rate-of-change will be set to zero
        i = np.array([self.dtype.y.names.index(k) for k in kwargs.keys()])
        v = np.array(kwargs.values())
        def clamped(t, y, ydot, f_data):
            """New RHS that prevents some elements from changing."""
            # Argh: clamped state variables may change slightly 
            # (within solver precision) despite setting ydot[i] = 0 below.
            # Therefore, we must also set the state variables to their clamped 
            # values on every invocation.
            
            # # Desperate debugging:
            # for k, v in kwargs.items():
            #     val = y[self.dtype.y.names.index(k)]
            #     if val == v:  print k, "should be", v, ", and it is"
            #     else:         print k, "should be", v, "but is", val
            
            y_ = np.copy(y)  # CVODE forbids modifying the state directly...
            y_[i] = v  # ...but may modify state variables even if ydot is always 0
            self.f_ode(t, y_, ydot, f_data)
            ydot[i] = 0
            return 0
        
        # Initialize clamped state variables.
        y = np.array(self.y).view(self.dtype.y)
        for k, v in kwargs.items():
            y[k] = v
        
        # Use original options when rerunning the Cvodeint initialization.
        oldkwargs = dict((k, getattr(self, k)) 
            for k in "chunksize maxsteps reltol abstol".split())
        
        clamped = Namedcvodeint(clamped, self.t, y, self.pr, **oldkwargs)
        
        # Disable any hard-coded stimulus protocol
        if "stim_amplitude" in clamped.dtype.p.names:
            clamped.pr.stim_amplitude = 0
        
        try:
            yield clamped # enter "with" block
        finally:
            for k in clamped.dtype.y.names:
                if k in self.dtype.y.names:
                    setattr(self.yr, k, getattr(clamped.yr, k))
    
    def rates(self, t, y, par=None):
        """
        Compute rates for a given state trajectory.
        
        Unfortunately, the CVODE machinery does not offer a way to return rates 
        during integration. This function re-computes the rates at each time 
        step for the given state.
        
        >>> vdp = Namedcvodeint()
        >>> t, y, flag = vdp.integrate()
        >>> ydot = vdp.rates(t, y)
        >>> ydot[0], ydot[-1]
        (rec.array([(0.0, 2.0)], dtype=[('x', '<f8'), ('y', '<f8')]), 
         rec.array([(0.780218..., 0.513757...)], dtype=...))
        """
        t = np.atleast_1d(t).astype(float)
        y = np.atleast_2d(y).view(float)
        imax = len(t)
        ydot = np.zeros_like(y)
        with self.autorestore(_p=par):
            for i in range(len(t)):
                self.f_ode(t[i], y[i], ydot[i], None)
        ydot = ydot.squeeze().view(self.dtype.y, np.recarray)
        return ydot

class Recarraylink(object):
    """
    Dynamic link between a Numpy recarray and any array-like object.
    
    Example: Use variable names to access elements of a Sundials state vector.
    We need to specify the Numpy dtype (data type), which in this case is 
    built by the :class:`~cellmlmodels.cellmlmodel.Cellmlmodel` constructor.
    
    >>> vdp = Namedcvodeint() # default example: van der Pol model
    >>> ral = Recarraylink(vdp.y, vdp.dtype.y)
    >>> ral
    Recarraylink([-2.0, 0.0], [('x', '<f8'), ('y', '<f8')])
    
    Now, a change made to either object is mirrored in the other:
    
    >>> ral.x = 123
    >>> vdp.y[1] = 321
    >>> vdp.y
    [123.0, 321.0]
    >>> ral
    Recarraylink([123.0, 321.0], [('x', '<f8'), ('y', '<f8')])
    
    Fixed bug: Previously, a direct modification of the array-like object was
    not applied to the recarray if the next Recarraylink operation was setting
    a named field. (The old value of the recarray was modified with the new
    field, and the direct modification was just forgotten. Afterwards, the
    wrong result overwrote the array-like object.)
    
    Here's an example relevant to Pysundials and CVODE.
    
    >>> y0 = np.array([0.0, 1.0, 2.0]) # default initial state
    >>> y = np.zeros(len(y0)) # state vector
    >>> ral = Recarraylink(y, dict(names=list("abc"), formats=[float]*len(y)))
    >>> y[:] = y0 # re-initialize state vector
    >>> ral.b = 42 # external modification
    >>> ral # reflects both changes
    Recarraylink([  0.  42.   2.], [('a', '<f8'), ('b', '<f8'), ('c', '<f8')])
    
    (Before the bugfix, output was
    Recarraylink([  0.  42.   0.], [('a', '<f8'), ('b', '<f8'), ('c', '<f8')])
    not reflecting the direct modification.)
    
    .. todo:: Recarraylink sometimes fails, perhaps for clamped models where I 
       change the NVector and cvode_mem.
    """
    def __init__(self, x, dtype):
        """Constructor. x : array-like object. dtype : numpy data type."""
        # use "object" class to avoid calling __getattr__ and __setattr__
        # before we're ready
        object.__setattr__(self, "_x", x)
        object.__setattr__(self, "_xa", np.array(x))
        object.__setattr__(self, "_xr", self._xa.view(dtype).view(np.recarray))
    
    def __getattr__(self, attr):
        """Return recarray field after copying values from array-like object"""
        try:
            self._xa[:] = self._x # from array-like object to ndarray view 
            return self._xr[attr] # access attribute via recarray view
        except ValueError: # pass through attributes not in self._xr
            super(Recarraylink, self).__getattr__(attr)
    
    def __setattr__(self, attr, val):
        """Set recarray field and copy values to array-like object"""
        try:
            self._xa[:] = self._x # refresh ndarray in case of extern chg to _x
            self._xr[attr] = val # set attribute of recarray view
            self._x[:] = self._xa # copy ndarray view's values to array-like obj
        except ValueError: # pass through attributes not in self._xr
            super(Recarraylink, self).__setattr__(attr, val)
    
    def __array__(self):
        """
        Underlying record array for array operations.
        
        >>> np.copy(Recarraylink([-2.0, 0.0], [('x', '<f8'), ('y', '<f8')]))
        array([(-2.0, 0.0)], 
              dtype=[('x', '<f8'), ('y', '<f8')])
        
        Without __array__, output would have been
        array(Recarraylink([-2.0, 0.0], [('x', '<f8'), ('y', '<f8')]), 
            dtype=object)
        """
        return self._xr
    
    def __repr__(self):
        """
        Detailed string representation.
        
        >>> Recarraylink([123.0, 321.0], [('x', '<f8'), ('y', '<f8')])
        Recarraylink([123.0, 321.0], [('x', '<f8'), ('y', '<f8')])
        """
        return "Recarraylink(%s, %s)" % (self._x, self._xr.dtype)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=
        doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS)
