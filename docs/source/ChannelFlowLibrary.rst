The Channel Flow Library
========================
The channel flow library contains Modelica building blocks to compose water system models from. The modeling objects are grouped into different categories. There are three main categories that would be used directly by the modeller: hydraulic, simple routing and salt. Hydraulic blocks model the open water channel flow with spatially discretized nodes. Simple routing has elements that has only one water level, or only volume. The salt block is an experimental version to model dispersive and convective salt transport.

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
This group of modeling objects contains building blocks to compose water balance models. Blocks have only one water level or only volume.

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

The connector ``QOut.Q`` is free (input), therefore it is typically connected to a time series which represents a as control variable. In optimization models, the control variable is an optimization variable, in simulation models the value of a control variable is computed according to the feedback control logic, which is typically specified in Python code. 

Steady
""""""

The *Steady* branch basically passes inflow to outflow. Lateral inflow and forcing can be added. This model object is used to model river sections where travel time and wave damping can be neglected. Typical application cases are reservoir models and water balance models. 

Nodes
^^^^^

Node
""""

A *Node* connects multiple branches. It can be used to model bifurcations and confluences in the water system model. It is designed such that the first inputs are the known ones and the last one is calculated from mass balance.

Example:

    ``Deltares.ChannelFlow.SimpleRouting.Nodes.Node Alder(nin = 2, nout = 1, n_QForcing = 0, QIn.Q(each nominal = 0.3), QOut.Q(each nominal = 10))``

In this example, the *Node*  with name *Alder* has two inflow connectors (``nin``) and one outflow connector (``nout`` ) to represent a confluence. This node is used in an optimization model and nominal values are assigned to the inflow and outflow connectors. It is possible to specify forcing outflow (``n_QForcing``) to a *Node*. Typically, forcing is used to model extractions from a channel network. 


Reservoir
^^^^^^^^^

The *Reservoir* node is used for modelling reservoirs. The basic equation is the storage equation

.. math:: \frac{\partial V}{\partial t} = Q_\mathrm{in} - Q_\mathrm{out} + Q_\mathrm{forcing} + Q_\mathrm{lateral}

The reservoir outflow is split into turbine flow and spill:

.. math:: Q_\mathrm{out} = Q_\mathrm{turbine} + Q_\mathrm{spill}

The basic idea is to distinguish between the total reservoir release that is available for hydropower production and the portion that is not used for hydropower production (spill). An optimization model usually has a goal to minimize the spill such that water preferably is guided through the turbine. Note that the spill flow component summarizes all flow components that are not turbine flow. Usually this is the flow through the bottom outlet and the spillway flow. 

:math:`Q_\mathrm{forcing}` can be used to represent other extractions from the reservoir, for example extractions for drinking water supply.  :math:`Q_\mathrm{lateral}` typically represents inflow into the reservoir from other sources than the main inflow. 

Storage
^^^^^^^

A *Storage* node represents a more general storage in a model. The basic equation is:

.. math:: \frac{\partial V}{\partial t} = Q_\mathrm{in} - Q_\mathrm{out} + Q_\mathrm{forcing}

While for a reservoir node the inflow is often given as external inflow time series or depends on upstream water management, for a *Storage* node often both the inflow `QIn.Q` and the outflow `QOut.Q` are control variables and either determined within the optimization (optimization model) or by a feedback control logic (simulation model). Thus the only difference between the reservoir and the storage blocks (except for the reservoir allowing for lateral flows) that the reservoir block has two free outflows: release and spill - thus both of them will be set by optimization or should be set by simulation. 

Structures
^^^^^^^^^^

DischargeControlledStructure
""""""""""""""""""""""""""""

The *DischargeControlledStructure* takes a discharge value. Typically, this modeling object is used to represent a hydraulic structure like a weir or a pump and apply the discharge to it. The value of the discharge `QOut.Q` is the control variable. 


Hydraulic
---------
This group of modeling objects contains building blocks for hydraulic modelling, basically one-dimensional open channel flow. 


Salt
----
This group contains building blocks for modeling salt water intrusion.


