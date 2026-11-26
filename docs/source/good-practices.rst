Good modelling practices
========================

The principal modelling practices identified from the example fixes are
summarized below. The general objective is to create Modelica models that are
mathematically well-posed, numerically robust, and valid according to the
Modelica Language Specification. A model should not depend on permissive or
tool-specific behaviour.

Avoid over- and under-specification
------------------------------------

Every equation should constrain an actual unknown.

* Do not initialize boundary variables that are already determined by boundary
  conditions or input equations.
* Do not introduce redundant initial equations.
* Initialize only genuine dynamic states and degrees of freedom.
* Ensure that the system contains enough independent equations for its
  unknowns.
* Treat solver warnings about insufficient degrees of freedom as modelling
  errors, even when the solver obtains a solution through a fallback.

For example, initialize only the internal channel flows when the boundary flows
are already prescribed:

.. code-block:: modelica

   Channel.Q[2:Channel.n_level_nodes] =
       fill(Inflow_Q, Channel.n_level_nodes - 1);

Initializing ``Q[1]`` and ``Q[n_level_nodes + 1]`` as well would repeat
information already provided by the boundary equations.

As a rule of thumb, initialize states rather than variables whose values
already follow algebraically from other equations.

Keep external inputs genuinely external
---------------------------------------

A top-level input supplied externally by RTC-Tools should not have a binding
equation in its declaration.

Instead of:

.. code-block:: modelica

   input Real Inflow_Q(fixed = true) = Discharge.Q;

use an unbound input and a separate model equation:

.. code-block:: modelica

   input Real Inflow_Q(fixed = true);

   equation
     Inflow_Q = Discharge.Q;

A binding on a top-level input changes its Modelica semantics. An unbound
top-level input represents a value supplied from outside the model, which is
the intended behaviour for RTC-Tools time-series inputs and control variables.

Provide start values for unbound parameters
-------------------------------------------

Parameters that are assigned by RTC-Tools at run time should still have
meaningful ``start`` values when required:

.. code-block:: modelica

   parameter Real theta(start = 1.0);
   parameter Modelica.Units.SI.Duration step_size(start = 0.0);

RTC-Tools can override these values during execution, while standalone
Modelica tools can still instantiate and check the model.

The start value should have a clear physical or mathematical meaning. For
example:

* ``theta = 1.0`` represents the full nonlinear model.
* ``step_size = 0.0`` represents the intended default formulation.

Respect Modelica types
----------------------

Do not rely on implicit conversions accepted by a particular compiler.

For example, do not bind an ``Integer`` parameter to a ``Real`` literal:

.. code-block:: modelica

   parameter Integer delay_hours = 24.0;

Use an integer literal instead:

.. code-block:: modelica

   parameter Integer delay_hours = 24;

``Integer`` and ``Real`` are distinct Modelica types.

Use resolvable and preferably fully qualified type names. Instead of relying on
an undefined ``SI`` alias:

.. code-block:: modelica

   input SI.Stress reach_1_stress_u;

use:

.. code-block:: modelica

   input Modelica.Units.SI.Stress reach_1_stress_u;

The shorter form is appropriate only if ``SI`` has been explicitly imported.

Match array dimensions explicitly
---------------------------------

Do not rely on implicit scalar-to-array broadcasting. Both sides of a Modelica
equation should have compatible dimensions.

Instead of:

.. code-block:: modelica

   der(Channel.Q) = 0;

construct an array of the correct size:

.. code-block:: modelica

   der(Channel.Q) =
       fill(0.0, Channel.n_level_nodes + 1);

Similarly, use:

.. code-block:: modelica

   downstream.Q[2:downstream.n_level_nodes + 1] =
       fill(Inflow_Q, downstream.n_level_nodes);

rather than assigning a scalar directly to an array slice.

Scale optimization goals appropriately
--------------------------------------

Objective functions and constraints should be scaled according to their actual
units and expected physical magnitudes.

For example, if ``Q_pump`` is measured in ``m3/s`` and RTC-Tools time is
expressed in seconds, then ``der(Q_pump)`` has units of ``m3/s2``. A nominal
change of ``5 m3/s`` per hour should therefore be specified as:

.. code-block:: python

   function_nominal = 5.0 / 3600.0  # m3/s2

Using ``5.0`` directly would represent ``5 m3/s2``, which is a substantially
different scale. Poor scaling can cause objective gradients to approach the
solver tolerance and can make results dependent on where the solver happens
to terminate.

Keep post-processing outside the physical model
-----------------------------------------------

Modelica should describe physical relationships at each point in time.
Operations across the complete optimization horizon generally belong in the
RTC-Tools Python layer.

For example, a total over all time steps should be calculated during
post-processing:

.. code-block:: python

   total_revenue = np.sum(results["TotalSystemRevenue"])

It should not be represented as an artificial Modelica variable using array
operations on a quantity that is scalar at each time instant.

A useful division of responsibilities is:

* **Modelica:** physical relationships and component behaviour at each time.
* **RTC-Tools/Python:** optimization goals, horizon-wide calculations,
  post-processing, aggregation, and reporting.

Only add meaningful symbolic constraints
----------------------------------------

Do not add an optimization constraint if its expression has already reduced to
a constant.

For example, ``der_at()`` can return a numerical value when the corresponding
quantity is fully determined by input data. Such an expression has nothing
left to constrain and should not be added as a CasADi constraint.

A constant equality such as ``0 = 0`` provides no information. A constant
inconsistent equality makes the problem infeasible without constraining a
decision variable.

Write portable Modelica
-----------------------

A model should conform to the Modelica Language Specification rather than only
being accepted by a particular compiler.

Common portability issues include:

* scalar-to-array broadcasting;
* assigning ``Real`` literals to ``Integer`` parameters;
* unqualified or unresolved type names;
* binding externally supplied top-level inputs;
* omitting required parameter start values;
* using invalid array operators; and
* using different identifiers at the beginning and end of a class.

For example:

.. code-block:: modelica

   model ExampleOptimization
     extends Example;
   end ExampleOptimization;

The identifier following ``end`` should match the class identifier.

Checking models with both pymoca and a second Modelica implementation, such as
OpenModelica, is therefore a useful conformance test.

Test structural and numerical correctness
-----------------------------------------

A successful solver termination is not sufficient evidence that the model is
correct.

Tests should detect structural solver messages such as:

.. code-block:: text

   Too few degrees of freedom

After a structural model correction, verify that the intended numerical
behaviour is preserved. Depending on the change, this can include comparing:

* the number of variables and constraints;
* constraint Jacobian ranks;
* variable and constraint bounds;
* solver iteration logs;
* objective values;
* exported time series;
* reference output files; and
* results from source-based and flattened model compilation.

Small numerical differences should be interpreted in relation to solver
tolerances, scaling, and convergence behaviour. Reference output should only
be updated after confirming that the difference is understood and expected.

Checklist
---------

Use the following checklist when creating or reviewing an RTC-Tools Modelica
example:

#. Keep the equation system structurally balanced.
#. Initialize only genuine dynamic states and degrees of freedom.
#. Do not initialize variables already fixed by equations or boundary
   conditions.
#. Keep externally supplied top-level inputs unbound.
#. Give externally assigned parameters meaningful ``start`` values.
#. Match ``Integer``, ``Real``, ``Boolean``, and other types exactly.
#. Match array dimensions explicitly and use ``fill()`` where necessary.
#. Use resolvable, preferably fully qualified, Modelica type names.
#. Scale goals and constraints according to their actual units and expected
   magnitudes.
#. Keep horizon-wide calculations and reporting in Python when they are not
   physical Modelica equations.
#. Do not add constant expressions as CasADi or RTC-Tools constraints.
#. Test models with both RTC-Tools/pymoca and a standards-oriented Modelica
   compiler such as OpenModelica.
#. Treat structural solver warnings and fallbacks as model problems.
#. Verify numerical equivalence after structural changes.