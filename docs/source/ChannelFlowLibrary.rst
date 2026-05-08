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

Hydraulic
---------
This group of modeling objects contains building blocks for hydraulic modelling, basically one-dimensional open channel flow. 


Salt
----
This group contains building blocks for modeling salt water intrusion.


