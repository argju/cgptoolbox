"""
Pysundials (CVODE) wrapper for Python code autogenerated from cellml.org

..  plot::
    
    from cgp.physmod.cellmlmodel import Cellmlmodel
    plt.title("Van der Pol Heart cell model")
    vdp = Cellmlmodel()
    t, y, flag = vdp.integrate(t=[0, 20])
    plt.plot(t, y.x, t, y.y)
"""
from __future__ import with_statement # allows garbage collection of files
from urllib import urlopen # to download from CellML Code Generation Service
from contextlib import contextmanager
from collections import namedtuple
import os # to remove old compiled version of model modules

from nose.tools import nottest

from ..cvodeint.namedcvodeint import Namedcvodeint
from ..utils.commands import getstatusoutput
from ..utils.dotdict import Dotdict
from ..cvodeint import Cvodeint # CVODE wrapper
import numpy as np
from numpy import recarray # recarray allows named columns: y.V etc.
from cgp import physmod as cellmlmodels
import sys
import warnings
from ..utils.ordereddict import OrderedDict

__all__ = ["Cellmlmodel"]

# Ensure that cellmlmodels/cellml2py/ is a valid package directory
# This makes it easy to force re-generation of code by renaming cellml2py/
try:
    import cellml2py  # pylint: disable=W0611
except ImportError:
    # deferred import to minimize dependencies
    from ..utils.write_if_not_exists import write_if_not_exists
    dirname, _ = os.path.split(__file__)
    with write_if_not_exists(os.path.join(dirname, "cellml2py", "__init__.py")):
        pass # just create an empty __init__.py file

Legend = namedtuple("Legend", "name component unit")

def parse_variable(s):
    """
    Parse the generated legend string for a CellML variable.
    
    >>> parse_variable("x in component Main (dimensionless)")
    Legend(name='x', component='Main', unit='dimensionless')
    """
    name, s1 = s.split(" in component ")
    component, s2 = s1.rsplit(" (")
    unit, s3 = s2.rsplit(")")
    return Legend(name, component, unit)

def parse_legend(legend):
    """
    Parse the entry for each CellML variable in a legend list.
    
    >>> parse_legend(["x in component A (u)", "y in component B (v)"])
    Legend(name=('x', 'y'), component=('A', 'B'), unit=('u', 'v'))
    """
    L = []
    for i, s in enumerate(legend):
        if s:
            L.append(parse_variable(s))
        else:
            # This will be used as a Numpy field name and PyTables column name. 
            # The latter cannot start with double underscore in PyTables < 2.2.
            # http://www.pytables.org/trac/ticket/291
            leg = Legend(name="empty__%s" % i, component="", unit="")
            warnings.warn("Empty legend entry: %s" % (leg,))
            L.append(leg)
    L = Legend(*zip(*L))
    return L

def legend(model):
    """
    Parse the legends of a CellML model.
    
    >>> vdp = Cellmlmodel()
    >>> legend(vdp.model) # doctest: +NORMALIZE_WHITESPACE
    OrderedDict([('y', Legend(name=('x', 'y'), component=('Main', 'Main'), 
                    unit=('dimensionless', 'dimensionless'))), 
                 ('a', None), 
                 ('p', Legend(name=('epsilon',), component=('Main',), 
                    unit=('dimensionless',)))])
    
    The legend is available as an attribute of a Cellmlmodel object.
    Zipping it can be convenient sometimes.
    
    >>> zip(*vdp.legend["y"])
    [('x', 'Main', 'dimensionless'), ('y', 'Main', 'dimensionless')]
    """
    # Use OrderedDict and named tuples rather than Numpy record arrays
    # until we know the length of all strings.
    states, algebraic, voi, constants = model.createLegends()
    legend = OrderedDict([("y", states), ("a", algebraic), ("p", constants)])
    L = [(k, parse_legend(v) if v else None) for k, v in legend.items()]
    return OrderedDict(L)

ftype = np.float64 # explicit type declaration, can be used with cython

def dup(L):
    """
    Return duplicated elements of a list.
    
    >>> dup(list("aba"))
    ['a']    
    """
    uniq = set(L)
    L = list(L[:]) # copy the list to avoid side effects
    for i in uniq:
        L.remove(i)
    return L

def dtypes(legend):
    """
    Return Numpy data types for a CellML model; object with attributes y, p, a.   
    
    The result is a Dotdict, whose keys can be used as attributes.

    >>> from cgp.physmod.cellmlmodel import Cellmlmodel
    >>> vdp = Cellmlmodel()
    >>> d = dtypes(legend(vdp.model))
    
    Now, d.y, d.p, d.a are dtype for state, parameters or "algebraic" variables.
    
    >>> d.y
    dtype([('x', '<f8'), ('y', '<f8')])
    
    >>> d
    Dotdict({'a': None,
     'p': dtype([('epsilon', '<f8')]),
     'y': dtype([('x', '<f8'), ('y', '<f8')])})
    """
    # Dict of duplicates in each legend item (states, algebraics, parameters)
    d = dict((k, dup(v.name)) for k, v in legend.items() if v)
    d = dict((k, v) for k, v in d.items() if v) # drop items with no duplicates
    if d: # duplicate names, must be disambiguated
        import copy
        L = copy.deepcopy(legend)
        warnings.warn("Duplicate names: %s" % d)
        for k, dupnames in d.items():
            name = list(L[k].name) # convert tuple to mutable list
            # append __i to each duplicate item, where i is its index
            for i, n in enumerate(name):
                if n in dupnames:
                    name[i] = "%s__%s" % (n, i)
            assert not dup(name), "Failed to disambiguate names: %s" % d
            L[k] = L[k]._replace(name=tuple(name))
    else: # no duplicates, use legend as is
        L = legend
    # Generate data type for each legend item
    DT = [(k, np.dtype([(n, ftype) for n in v.name]) if v else None) 
        for k, v in L.items()]
    return Dotdict(DT)

#: Python code appended to that which is autogenerated from CellML
py_addendum = '''
### Added by cellmlmodel.py ###

# @todo: The following module-level variables are shared across instances.
#        It might be better to wrap them in a class, allowing each instance of 
#        the same model to have its own parameter vector.

import sys
import numpy as np

ftype = np.float64 # explicit type declaration, can be used with cython
y0 = np.zeros(sizeStates, dtype=ftype)
ydot = np.zeros(sizeStates, dtype=ftype)
p = np.zeros(sizeConstants, dtype=ftype)
algebraic = np.zeros(sizeAlgebraic, dtype=ftype)

y0[:], p[:] = initConsts()

exc_info = None

# Sundials calling convention: https://computation.llnl.gov/casc/sundials/documentation/cv_guide/node6.html#SECTION00661000000000000000

def ode(t, y, ydot, f_data):
    """
    Compute rates of change for differential equation model.
    
    Rates are written into ydot[:]. 
    f_data is ignored, but required by the CVODE interface.
    
    The function returns 0 on success and -1 on failure.
    
    >>> ode(None, None, None, None)
    -1
    
    For debugging in case of failure, exception info is stored in the 
    module-level variable exc_info. (The message ends in "unsubscriptable" 
    under Python 2.6 but "not subscriptable" under Python 2.7, hence the 
    ellipsis.) Unfortunately, this is currently not implemented in a compiled 
    ODE. It will check the type of arguments before executing, but I am not 
    sure what happens in case of run-time errors inside the ODE.
    
    >>> exc_info
    (<type 'exceptions.TypeError'>,
    TypeError("'NoneType' object is ...subscriptable",),
    <traceback object at 0x...>)
    """
    global exc_info
    exc_info = None
    try:
        ydot[:] = computeRates(t, y, p)
        return 0
    except StandardError:
        import sys
        exc_info = sys.exc_info()
        return -1

def rates_and_algebraic(t, y):
    """
    Compute rates and algebraic variables for a given state trajectory.
    
    Unfortunately, the CVODE machinery does not offer a way to return rates and 
    algebraic variables during integration. This function re-computes the rates 
    and algebraics at each time step for the given state.
    
    This returns a simple float array; 
    :meth:`cgp.physmod.cellmlmodel.Cellmlmodel.rates_and_algebraic`
    will cast them to structured arrays with named fields.
    
    This version is pure Python; 
    :func:`~cgp.physmod.cythonize.cythonize`
    will generate a faster version.
    """
    imax = len(t)
    # y can be NVector, unstructured or structured Numpy array.
    # If y is NVector, its data will get copied into a Numpy array.
    y = np.array(y).view(float)
    ydot = np.zeros_like(y)
    alg = np.zeros((imax, len(algebraic)))
    for i in range(imax):
        ydot[i] = computeRates(t[i], y[i], p)
        if len(algebraic):
            # need np.atleast_1d() because computeAlgebraic() uses len(t)
            alg[i] = computeAlgebraic(p, y[i], np.atleast_1d(t[i])).squeeze()
    return ydot, alg
'''


class Cellmlmodel(Namedcvodeint):
    """
    Class to solve CellML model equations.

    ..  plot::
    
        vdp = Cellmlmodel() # default van der Pol model
        t, Yr, flag = vdp.integrate(t=[0, 20])
        plt.plot(t,Yr.view(float))
    
    The constructor will download autogenerated Python code from cellml.org if
    possible, and otherwise look for a corresponding .py.orig file in the
    cgp/physmod/cellml2py/ directory. This code is wrapped to be compatible with
    CVode, and saved as a .py file in the same directory.
    
    If ``use_cython=True`` (the default), the code is rewrapped for `Cython
    <http://www.cython.org>`_ and compiled for speed. Compiled models reside in
    cgp/physmod/cellml2py/cython/ where each compiled model has a subdirectory
    containing several files. The .pyx file and setup.py file can be tweaked by
    hand if required, and manually recompiled by changing to that directory and
    running python setup.py build_ext --inplace The compiled module has
    extension .so (Linux) or .pyd (Windows).
        
    Compiling the model causes some minor differences in behaviour, see
    :func:`~cgp.test.test_cellmlmodel.test_compiled_behaviour` for details.
    """
    
    def __init__(self,  # pylint: disable=W0102
        exposure_workspace="5756af26cfb20a7f66a51f66af10a70a/vanderpol_vandermark_1928",
        urlpattern="http://models.cellml.org/exposure/%(exposure)s/%(workspace)s.cellml/@@cellml_codegen/Python/raw", 
        t=[0, 1], y=None, p=None, purge=False, rename={}, use_cython=True, 
        **kwargs):
        """
        Wrap autogenerated CellML->Python for use with pysundials

        M = Cellmlmodel(exposure_workspace) downloads/caches Python code
        autogenerated for the CellML model identified by exposure/workspace,
        and wraps it in a class with convenience attributes and functions for
        use with pysundials.

        If the non-wrapped Python code is in a local file, 
        e.g. exported from OpenCell, http://www.cellml.org/tools/opencell/ 
        use the "file://" protocol.
        
        >>> newmodel = Cellmlmodel("/newmodel", "file:c:/temp/exported.py")
        ... # doctest: +SKIP
        
        Here, "newmodel" is whatever name you'd like for the wrapper module, 
        and "exported.py" is whatever name you saved the exported code under.
        (Strictly speaking, the URL should be "file:///c:/temp/exported.py", 
        but the simpler version is also accepted by urllib.urlopen().)
        
        The constructor arguments are as follows:
        
        exposure_workspace: identifiers in the repository at cellml.org,
        e.g. "732c32162c845016250f234416415bfc7601f41c/vanderpol_vandermark_1928_version01"
        for http://models.cellml.org/exposure/2224a49c6b39087dad8682648658775d
        
        urlpattern : URL to non-wrapped Python code for model, with
        %(workspace)s and %(exposure)s placeholders for e.g.
        732c32162c845016250f234416415bfc7601f41c
        vanderpol_vandermark_1928_version01
        
        t, y : as for Cvodeint
        
        p : optional parameter vector
        
        purge : (re-)download model even if the file is already present?
        
        rename : e.g. dict with possible keys "y", "p", "a", whose values are 
        mappings for renaming variables. You should rarely need this, but it is 
        useful to standardize names of parameters to be manipulated, see e.g. 
        ap_cvode.Tentusscher.__init__().
        
        use_cython: if True, wrap the model for Cython and compile.
        Cython files are placed in cgp/physmod/cellml2py/cython/modulename/, 
        and cgp.physmod.cellml2py.cython.modulename.modulename is used 
        in place of cgp.physmod.cellml2py.modulename.
        
        >>> Cellmlmodel().dtype
        Dotdict({'a': None,
         'p': dtype([('epsilon', '<f8')]),
         'y': dtype([('x', '<f8'), ('y', '<f8')])})
        >>> Cellmlmodel(rename={"y": {"x": "V"}, "p": {"epsilon": "newname"}}).dtype
        Dotdict({'a': None,
         'p': dtype([('newname', '<f8')]),
         'y': dtype([('V', '<f8'), ('y', '<f8')])})
        
        See class docstring: ?Cellmlmodel for details.
        """
        exposure, workspace = exposure_workspace.split("/", 1)
        self.name = modelname = workspace + exposure
        modelfilename = cellmlmodels.__path__[0]
        modelfilename += '/cellml2py/' + modelname + ".py"
        modulename = 'cgp.physmod.cellml2py.' + modelname
        # import the model module, creating it if necessary
        try:
            if purge:
                try:
                    os.remove(modelfilename + "c")
                except OSError:
                    pass
                raise ImportError
            __import__(modulename)
            try:
                with open(modelfilename, "rU") as f:
                    self.py_code = f.read()
            except IOError:
                self.py_code = "Source file open failed"
            try:
                with open(modelfilename + ".orig", "rU") as f:
                    self.py_code_orig = f.read()
            except IOError:
                self.py_code_orig = "Source file open failed"

        except ImportError:
            # try to download Python code autogenerated from CellML
            url = urlpattern % locals()
            origfile = modelfilename + ".orig"
            try:
                py_code = "".join(urlopen(url).readlines())
            except IOError:
                py_code = ""
            if "computeRates" in py_code:
                # save original file for later reference, e.g. cythonization
                with open(origfile, "w") as f:
                    f.write(py_code)
            else:
                # look for existing .py.orig file in cgp/physmod/cellml2py
                msg = ("Failed to get autogenerated Python code. " +
                    "URL did not provide valid Python code. %s\n" + 
                    "If you have a .cellml file, try exporting it to Python " +
                    "using OpenCell and save it under the 'local file' name " +
                    "given below.\nOpenCell: " +
                    "http://www.cellml.org/tools/downloads/opencell/releases/" +
                    "\nURL: " + url + "\nLocal file: " + origfile)
                try:
                    with open(origfile) as f:
                        py_code = f.read()
                except IOError:
                    raise IOError(msg % "Local file does not exist.")
                assert "computeRates" in py_code, (msg % 
                    "Local file does exist, but Python code is not valid.")
            
            self.py_code = self.py_code_orig = py_code
            # prevent style checking of autogenerated code
            self.py_code = "# pylint: disable-all\n" + self.py_code
            # ensure floating-point division even of integer-valued parameters
            self.py_code = "from __future__ import division\n" + self.py_code
            # add ode(y, t, p) wrapper and a few lines of initialization
            self.py_code += py_addendum
            self.py_code += "\n# Downloaded from %s\n" % url
            # use numpy.core.* instead of math.* to work with arrays
            self.py_code = self.py_code.replace("from math import *", 
                "from numpy.core import *")
            # fix bug in CellML code generation
            self.py_code = self.py_code.replace("VOI", "voi")
            # write generated code to a module file
            with open(modelfilename, "w") as f:
                f.write(self.py_code)
            __import__(modulename)
        
        # store a reference to the model module
        if use_cython:
            self.model = self.cythonize(modelname, modulename, modelfilename)
        else:
            self.model = sys.modules[modulename]
        
        if y is None:
            y = self.model.y0
        self.legend = legend(self.model)
        dtype = dtypes(self.legend)
        # Rename fields if requested
        for i in "a", "y", "p":
            if i in rename:
                L = eval(str(dtype[i]))
                for j, (nam, typ) in enumerate(L):
                    if nam in rename[i]:
                        L[j] = (rename[i][nam], typ)
                dtype[i] = np.dtype(L)
        
        # if there are no parameters or algebraic variables, make empty recarray
        try:
            pr = self.model.p.view(dtype.p).view(recarray)
        except TypeError:
            pr = np.array([]).view(recarray)
        try:
            self.algebraic = self.model.algebraic.view(dtype.a).view(recarray)
        except TypeError:
            self.algebraic = np.array([]).view(recarray)
        self.y0r = self.model.y0.view(dtype.y).view(recarray)
        super(Cellmlmodel, self).__init__(self.model.ode, t, 
            y.view(dtype.y), pr, **kwargs)
        assert all(dtype[k] == self.dtype[k] for k in self.dtype)
        self.dtype.update(dtype)
        self.originals["y0r"] = self.y0r
        if p:
            self.model.p[:] = p
    
    @contextmanager
    def dynclamp(self, setpoint, R=0.02, V="V", ion="Ki", scale=None):
        """
        Derived model with state value dynamically clamped to a set point.
        
        Input arguments:
        
        * setpoint : target value for state variable
        * R=0.02 : resistance of clamping current
        * V="V" : name of clamped variable
        * ion="Ki" : name of state variable carrying the clamping current
        * scale=None : "ion" per "V", 
          default :math:`Acap * Cm / (Vmyo * F)`, see below
        
        Clamping is implemented as a dynamically applied current that is 
        proportional to the deviation from the set point::
        
            dV/dt = -(i_K1 + ... + I_app)
            I_app = (V - setpoint) / R
        
        Thus, if V > setpoint, then I_app > 0 and serves to decrease dV/dt.
        To account for the charge added to the system, a current proportional 
        to I_app is added to a specified ion concentration, by default Ki. This 
        needs to be scaled according to conductance and other constants.
        The default is as for the Bondarenko model::
         
            scale = Acap * Cm / (Vmyo * F)
            dKi/dt = -(i_K1 + ... + I_app) * scale
        
        Example with voltage clamping of Bondarenko model.
        
        >>> bond = Cellmlmodel("11df840d0150d34c9716cd4cbdd164c8/"
        ...     "bondarenko_szigeti_bett_kim_rasmusson_2004_apical")
        >>> with bond.dynclamp(-140) as clamped:
        ...     t, y, flag = clamped.integrate(t=[0, 10])
        
        The clamped model has its own instance of the CVODE integrator and 
        state NVector.
        
        >>> clamped.cvode_mem is bond.cvode_mem
        False
        
        However, changes in state are copied to the original on exiting the 
        'with' block.
        
        >>> "%7.2f" % bond.yr.V
        '-139.68'
        
        Unlike .clamp(), .dynclamp() does not allow you to change the setpoint 
        inside the "with" block. Instead, just start a new "with" block.
        (Changes to state variables remain on exit from the with block.)
        
        >>> with bond.dynclamp(-30) as clamped:
        ...     t0, y0, flag0 = clamped.integrate(t=[0, 10])
        >>> with bond.dynclamp(-10) as clamped:
        ...     t1, y1, flag1 = clamped.integrate(t=[10, 20])
        >>> np.concatenate([y0[0], y0[-1], y1[0], y1[-1]])["V"].round(2)
        array([-139.68,  -29.96,  -29.96,  -10.34])
        
        Naive clamping with dV/dt = 0 and unspecified clamping current, 
        like .clamp() does, equals the limit as R -> 0 and scale = 0.
        
        >>> with bond.autorestore(): # FIXME: bond.autorestore(V=0) doesn't work
        ...     bond.y[0] = 0 # FIXME: bond.yr.V = 0 does not work
        ...     with bond.dynclamp(-140, 1e-10, scale=0) as clamped:
        ...         t, y, flag = clamped.integrate(t=[0, 0.1])
        
        Although V starts at 0, it gets clamped to the setpoint very quickly.
        
        >>> y.V[0]
        array([   0.])
        >>> t[y.V.squeeze() < -139][0]
        4.94...e-10
        """
        if scale is None:
            p = self.pr
            scale = p.Acap * p.Cm / (p.Vmyo * p.F)
        
        # Indices to state variables whose rate-of-change will be modified.
        iV = self.dtype.y.names.index(V)
        iion = self.dtype.y.names.index(ion)
        
        def dynclamped(t, y, ydot, f_data):
            """New RHS that prevents some elements from changing."""
            self.f_ode(t, y, ydot, f_data)
            I_app = (y[iV] - setpoint) / R
            ydot[iV] -= I_app
            ydot[iion] -= I_app * scale
            return 0
        
        y = np.array(self.y).view(self.dtype.y)
        
        # Use original options when rerunning the Cvodeint initialization.
        oldkwargs = dict((k, getattr(self, k)) 
            for k in "chunksize maxsteps reltol abstol".split())
        
        clamped = Namedcvodeint(dynclamped, self.t, y, self.pr, **oldkwargs)
        
        # Disable any hard-coded stimulus protocol
        if "stim_amplitude" in clamped.dtype.p.names:
            clamped.pr.stim_amplitude = 0
        
        try:
            yield clamped # enter "with" block
        finally:
            for k in clamped.dtype.y.names:
                if k in self.dtype.y.names:
                    setattr(self.yr, k, getattr(clamped.yr, k))
    
    def rates_and_algebraic(self, t, y, par=None):
        """
        Compute rates and algebraic variables for a given state trajectory.
        
        Unfortunately, the CVODE machinery does not offer a way to return rates and 
        algebraic variables during integration. This function re-computes the rates 
        and algebraics at each time step for the given state.
        
        ..  plot::
        
            exposure_workspace = "b0b1820b1376263e16c6086ca64d513e/bondarenko_szigeti_bett_kim_rasmusson_2004_apical"
            bond = Cellmlmodel(exposure_workspace, t=[0, 20])
            bond.yr.V = 100 # simulate stimulus
            t, y, flag = bond.integrate()
            ydot, alg = bond.rates_and_algebraic(t, y)
            plt.plot(t, alg.J_xfer, '.-', t, y.Cai, '.-')
        """
        m = self.model
        t = np.atleast_1d(t).astype(float)
        y = np.atleast_2d(y)
        imax = len(t)
        # y = y.view(ftype) # done already in rates_and_algebraic
        with self.autorestore(_p=par):
            ydot, alg = m.rates_and_algebraic(t, y)
        # ydot = ydot.view(self.dtype.y, np.recarray).squeeze()
        # alg = alg.view(self.dtype.a, np.recarray).squeeze()
        ydot = ydot.squeeze().view(self.dtype.y, np.recarray)
        alg = alg.squeeze().view(self.dtype.a, np.recarray)
        return ydot, alg
    
    def cythonize(self, modelname, modulename, modelfilename):
        """
        Return Cython code for this model (further hand-tweaking may be needed).
        
        This just imports and calls 
        :func:`cgp.physmod.cythonize.cythonize_model`.
        """
        modulename_cython = modulename.replace(".%s" % modelname, 
            ".cython.%s.%s" % (modelname, modelname))
        try:
            __import__(modulename_cython)
            return sys.modules[modulename_cython]
        except ImportError:
            # deferred import to minimize dependencies
            from cgp.physmod.cythonize import cythonize_model
            from ..utils.write_if_not_exists import write_if_not_exists
            pyx, setup = cythonize_model(self.py_code_orig, modelname)
            pyxname = modelfilename.replace("%s.py" % modelname, 
                "cython/%s/%s.pyx" % (modelname, modelname))
            dirname, _ = os.path.split(pyxname)
            setupname = os.path.join(dirname, "setup.py")
            # make the cython and model directories "packages"
            cyinitname = os.path.join(dirname, os.pardir, "__init__.py")
            modelinitname = os.path.join(dirname, "__init__.py")
            with write_if_not_exists(cyinitname):
                pass # just create an empty __init__.py file
            with write_if_not_exists(modelinitname):
                pass # just create an empty __init__.py file
            with write_if_not_exists(pyxname) as f:
                f.write(pyx)
            with write_if_not_exists(setupname) as f:
                f.write(setup)
            cmd = "python setup.py build_ext --inplace"
            status, output = getstatusoutput(cmd, cwd=dirname)
            # Apparently, errors fail to cause status != 0.
            # However, output does include any error messages.
            if "cannot find -lsundials_cvode" in output:
                raise OSError("Cython-compilation of ODE right-hand side "
                    "failed because SUNDIALS was not found.\n"
                    "Status code: %s\nCommand: %s\n"
                    "Output (including errors):\n%s" % (status, cmd, output))
            if status != 0:
                raise RuntimeError("'%s'\nreturned status %s:\n%s" % 
                    (cmd, status, output))
            try:
                __import__(modulename_cython)
                return sys.modules[modulename_cython]
            except StandardError, exc:
                raise ImportError("Exception raised: %s: %s\n\n"
                    "Cython compilation may have failed. "
                    "The compilation command was:\n%s\n\n"
                    "The output of the compilation command was:\n%s"
                    % (exc.__class__.__name__, exc, cmd, output))
    
    @nottest
    def makebench(self):
        """
        Return IPython code for benchmarking compiled vs. uncompiled ode.
        
        >>> print Cellmlmodel().makebench()
        # Run the following from the IPython prompt:
        import os
        import cgp.physmod... as m
        reload(m)
        from utils.capture_output import Capture
        with Capture() as cap:
            print "##### Timing the current version #####"
            timeit m.ode(0, m.y0, m.ydot, None)
        ...
        
        Executing this in IPython gives something like this for a model whose 
        right-hand-side module has been compiled.
        ##### Timing the current version #####
        10000 loops, best of 3: 22.7 us per loop
        ##### Timing pure Python version #####
        100 loops, best of 3: 2.28 ms per loop        
        """
        template = """# Run the following from the IPython prompt:
import os
import %s as m
reload(m)
from utils.capture_output import Capture
with Capture() as cap:
    print "##### Timing the current version #####" 
    timeit m.ode(0, m.y0, m.ydot, None)

src = m.__file__
_, ext = os.path.splitext(src)
if ext in [".so", ".pyd"]:
    backup = src + ".backup"
    os.rename(src, backup)
    try:
        reload(m)
        with cap:
            print "##### Timing pure Python version #####"
            timeit m.ode(0, m.y0, m.ydot, None)
    
    finally:
        os.rename(backup, src)
        print "##### Restored compiled version #####"
        reload(m)

print cap
"""
        return template % self.model.__name__
    
if __name__ == "__main__":
    import doctest
    failure_count, test_count = doctest.testmod(optionflags=
        doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS | 
        doctest.REPORT_ONLY_FIRST_FAILURE)
    print """
        NOTE: You may see AttributeError when pysundials tries to __del__
        NVector objects that are None. This is probably not a problem.
        
        Also, bugs in the CellML code generation cause a few warnings for the 
        Bondarenko model. 
        """
    if failure_count == 0:
        print """
            All doctests passed.
            """
