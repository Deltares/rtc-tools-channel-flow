The Channel Flow Library
========================
The channel flow library contains Modelica building blocks to compose water system models from. The modeling objects are grouped into different categories. 

Media
-----
The default medium is fresh water, salt water with variable density is also supported for some model objects. 

Interfaces
----------
Interfaces define parameters and flow direction of an object. 

Internal
--------
Building blocks of partial models. These internal models make use of interfaces to define the connection properties of an object. 

SimpleRouting
-------------
This group of modeling objects contains building blocks to compose water balance models. 

BoundaryConditions
^^^^^^^^^^^^^^^^^^
Inflow
""""""
The *Inflow* object is used to set an inflow boundary condition. It expects a volume flow rate (discharge) assigned to the connector ``Q``. 

Example: ``TroutLake_Inflow.Q = 23.0`` sets a constant inflow to the *Inflow* object *TroutLake_Inflow* of 23.0. Instead of a number, an input time series can be set. 

Terminal
""""""""
The *Terminal* node is used as downstream boundary condition. 

Example: ``RiverCity.Q`` collects the discharge that flows out of the model domain at the node *RiverCity*. 

Branches
^^^^^^^^

Delay
"""""
A *Delay* branch delays the outflow with respect to the inflow in seconds. The duration should be a multiple of the time step to ensure a closed water balance. 

Integrator
""""""""""

The *Integrator* branch contains a storage element. The governing equation is

.. math:: \frac{\partial V}{\partial t} = Q_\mathrm{in} - Q_\mathrm{out} + Q_\mathrm{forcing} + Q_\mathrm{lateral}

The connector ``QOut.Q`` is typically connected to a time series which represents a as control variable. In optimization models, the control variable is an optimization variable, in simulation models the value of a control variable is computed according to the feedback control logic, which is typically specified in Python code. 

Steady
""""""

The *Steady* branch basically passes inflow to outflow, but lateral inflow and forcing can be added. 

Nodes
^^^^^

Node
""""

A *Node* connects multiple branches. It can be used to model bifurcations and confluences in the water system model. 

Example:

    ``Deltares.ChannelFlow.SimpleRouting.Nodes.Node Alder(nin = 2, nout = 1, n_QForcing = 0, QIn.Q(each nominal = 0.3), QOut.Q(each nominal = 10))``

In this example, the *Node*  with name *Alder* has two inflow connectors (``nin``) and one outflow connector (``nout`` ) to represent a confluence. This node is used in an optimization model and nominal values are assigned to the inflow and outflow connectors. It is possible to specify forcing outflow (``n_QForcing``) to a *Node*. Typically, forcing is used to model extractions from a channel network. 


Reservoir
^^^^^^^^^

Storage
^^^^^^^

Structures
^^^^^^^^^^

DischargeControlledStructure
""""""""""""""""""""""""""""

The *DischargeControlledStructure* takes a discharge value. Typically, this modeling object is used to represent a hydraulic structure like a weir or a pump and apply the discharge to it. The value of the discharge 


Hydraulic
---------
This group of modeling objects contains building blocks for hydraulic modelling, basically one-dimensional open channel flow. 


Salt
----
This group contains building blocks for modeling salt water intrusion.


