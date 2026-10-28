Node
====

A *Node* connects multiple branches. It can be used to model
  bifurcations and confluences in the water system model.

Example:

::

   Deltares.ChannelFlow.SimpleRouting.Nodes.Node Alder(
       nin = 2,
       nout = 1,
       n_QForcing = 0,
       QIn.Q(each nominal = 0.3),
       QOut.Q(each nominal = 10))

In this example, the *Node* with name *Alder* has two inflow
  connectors and one outflow connector (``nout``) to represent a confluence. This node is used in an optimization model and nominal values are assigned to the inflow and outflow connectors. It is possible to specify forcing outflow (``n_QForcing``) to a *Node*. Typically, forcing is used to model extractions from a channel network.
