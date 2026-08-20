import re

import numpy as np


class _ObjectParameterWrapper(object):
    """
    Python wrapper class for Modelica models/classes.

    Non-nested parameters in the model can be accessed as an attributes of the
    Python object. This is not the case for non-parameters, e.g. control
    variables or input time series.
    """

    def __init__(self, optimization_problem):
        self.optimization_problem = optimization_problem

    def _parse_array(self, optimization_problem, ks):
        # Figure out dimension of array
        inds = re.search(r"\[(.*?)\]$", ks[0])
        inds = re.findall(r"(\d+)", inds.group(1))
        n_dim = len(inds)

        pattern = r".*?\[" + ",".join([r"(\d+)"]*n_dim) + r"\]$"
        prog = re.compile(pattern)
        indices = [prog.match(x).groups() for x in ks]

        shape = np.zeros(n_dim, dtype=int)
        for i in range(n_dim):
            shape[i] = max((int(x[i]) for x in indices))

        arr = np.zeros(shape)

        parameters = optimization_problem.parameters(0)

        # TODO: Why are parameters stored as individual elements? Would be much easier to just to
        # parameters[k].getMatrixValue.toArray() or .toMatrix, to avoid looping/regex string parsing.
        for k in ks:
            inds = tuple(int(x)-1 for x in prog.match(k).groups())
            arr[inds] = float(parameters[k])

        return arr

    def __getattr__(self, attr):
        """
        Not found in regular class member variables or functions. Lookup in
        Modelica model's parameters.
        """
        try:
            # Array parameters are not stored as arrays, but as individual
            # elements. So if the parameter is an array, we will have to put
            # it back together again.
            ks = [x for x in self.optimization_problem.parameters(0).keys()
                  if x.startswith(self.symbol + "." + attr + "[")]
            if ks:
                return self._parse_array(self.optimization_problem, ks)
            else:
                return self.optimization_problem.parameters(0)[self.symbol + "." + attr]
        except KeyError:
            raise AttributeError
