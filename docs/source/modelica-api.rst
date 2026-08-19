Modelica API
============

Pumping Station
---------------

.. cpp:namespace:: Deltares::ChannelFlow::Hydraulic::Structures::PumpingStation

These Modelica blocks are part of the standard 
`ChannelFlow library <https://gitlab.com/deltares/rtc-tools-channel-
flow.git>`_. It consists of the following components:

Pump
  A pump model with a QHP (discharge, head, power) relationship, to be used
  for optimization of e.g. costs. It extends
  ``Deltares.ChannelFlow.Hydraulic.Structures.Pump``

Resistance
  Quadratic resistance.

PumpingStation
  Encapsulating class for Pump and Resistance objects.

.. note::

  ``Pump`` and ``Resistance`` objects should always be placed inside a
  ``PumpingStation`` object.

.. warning::

  Only the ``Pump`` and ``Resistance`` objects from this library should be placed inside a
  ``PumpingStation`` object. This means that the ``Pump`` class from
  ``Deltares.ChannelFlow.Hydraulic.Structures`` should not be used within ``PumpingStation``.
  It is allowed to use other ``ChannelFlow`` objects (e.g. ``Pump``, etc.) as upstream or 
  downstream connections to the ``PumpingStation``.

Pump
~~~~

.. cpp:class:: Pump : Deltares::ChannelFlow::Hydraulic::Structures::Pump

  Represents a single pump object. Because the power of the pump is seldom a
  linear function of `Q` and `H`, this class is wrapped by the Python API's
  :py:class:`~rtctools_hydraulic_structures.pumping_station_mixin.Pump` which
  turns the power equation specified by :cpp:var:`power_coefficients` into
  a set of inequality constraints:

  .. math::

    \begin{aligned} \
    P &\ge P_{min} \cdot S\\[5pt]
    P &\le P_{max} \cdot S\\[5pt]
    P &\ge C_{1,1} + C_{1,2} Q + \ldots - P_{max} \cdot \left(1 - S\right)\\[5pt]
    \end{aligned}

  where :math:`S` is the pump status (off = 0, on = 1), and :math:`P_{min}`
  and :math:`P_{max}` are the minimum and maximum possible power when the pump
  is on.

  With minimization of pumping costs (i.e. power), the optimization results
  will satisfy the first inequality constraint as if it was an equality
  constraint.

  .. cpp:var:: Real power_coefficients[_, _, _]

    The power coefficients describe the relationship between the discharge, head
    and power. For example, one can consider a fit of the pump power of the
    general form:

    .. math::

      P = C_{1,1} + C_{1,2} Q + C_{2,1} H + C_{2,2} Q H + C_{1,3} Q^2 + \ldots

    The power coefficients matrix corresponds to the coefficients `C` in the
    equation above. To guarantee that optimization finds a good and stable
    solution, we require the coefficients of this polynomial to be chosen such that
    the polynomial is convex over the entire domain.

    When several power coefficients matrices are given, the pump power will be
    greater than or equal to all resulting equations. In other words, the pump
    power will take the maximum value of all equations. The specification of
    multiple power coefficients can, for example, be used to define a fully
    linear power relation that can be handled by MILP solvers.

    .. note::

      Strictly speaking it would only have to be convex over the (automatically)
      extended working area domain, the size of which is not always known before
      run-time.

  .. cpp:var:: Real working_area[_, _, _]

    The working area array describes the polynomials bounding the convex set of
    allowed Q-H coordinates. These polynomials typically arise from one of the
    following properties:

    * Q-H curve at minimum pump speed
    * Q-H curve at maximum pump speed
    * Minimum required efficiency (e.g. 50%)
    * Minimum and/or maximum input power constraint
    * Cavitation constraints
    * NPSH constraints

    The first coordinate of the array is the polynomial number. For example,
    ``working_area[1, :, :]`` would describe the first working area polynomial.
    The order of Q and H coefficients is the same as in
    :cpp:var:`power_coefficients`.

  .. cpp:var:: Real working_area_direction[_]

    The polynomials in :cpp:var:`working_area` describe the polynomials, but do not
    yet indicate what side of this polynomial the Q-H combination has to be on. So
    for each of the polynomials in the working area we have to specify whether the
    expression should evaluate to a positive expression (`=1`), or a negative
    expression (`=-1`).

    .. note::

      It may become unnecessary to specify this in the future, if it is possible
      to figure out a way to determine this automatically based on the polynomials
      and their crossing points.

  .. cpp:var:: Integer head_option = 0

    What head to use for the pump head. This can be one of the following four options:

    \-1
      The upstream head

    0
      The differential head (i.e. downstream head minus upstream head)

    1
      The downstream head.

    2
      The differential head (i.e. downstream head minus upstream head) where the upsteam head is a determined by 
      the variable HW which can be set using a function of the users choice.

    .. note::

      If ``head_option=2`` is used, we assume that all pumps within the same pumping station
      also use ``head_option=2``. When using this head option the user must also provide a function
      that determines the value of HW at each time step.

  .. cpp:var:: Modelica::SIunits::Duration minimum_on = 0.0*3600

    The minimum amount of time in seconds a pump needs to be on before allowed
    to switch off again. This applies to all pumps in this pumping station.

    .. note::

      Only multiples of the (equidistant) time step are allowed.

  .. cpp:var:: Modelica::SIunits::Duration minimum_off = 0.0*3600

    The minimum amount of time in seconds a pump needs to be off before allowed
    to switch on again. This applies to all pumps in this pumping station.

    .. note::

      Only multiples of the (equidistant) time step are allowed.

  .. cpp:var:: Modelica::SIunits::Energy start_up_energy = 0.0

    The energy needed to start a pump. This will be multiplied with the energy
    price associated to each pump to calculate the costs.

  .. cpp:var:: Real start_up_cost = 0.0

    Costs in e.g. EUR or kg CO2 associated with a pump start up. Many pump
    switches could for example mean the pump life is shortened, or that more
    maintenance is needed. These could then be expressed in monetary value, and
    associated with pump start up.

    .. important::

      Make sure that the units of this value are of the same units as
      :cpp:var:`start_up_energy` times the energy price.

  .. cpp:var:: Modelica::SIunits::Energy shut_down_energy = 0.0

    Energy needed to shut down a pump. See equivalent parameter for pump start :cpp:var:`start_up_energy` for more information.

  .. cpp:var:: Real shut_down_cost = 0.0

    Cost associated with a pump shutdown. See equivalent parameter for pump start :cpp:var:`start_up_cost` for more information.


Resistance
~~~~~~~~~~

.. cpp:class:: Resistance : Deltares::ChannelFlow::Internal::HQTwoPort

  Represents a single quadratic resistance object relating the head loss to
  the discharge:

  .. math::

    dH = C \cdot Q^2

  Because a non-linear equality constraint is not allowed in convex
  optimization, this class is wrapped by the Python API's
  :py:class:`~rtctools_hydraulic_structures.pumping_station_mixin.Resistance`
  which turns it into two inequality constraints:

  .. math::

    \begin{aligned} \
    dH &\ge C \cdot Q^2 \\[5pt]
    dH &\le \alpha \cdot Q \
    \end{aligned}

  where the second constraint makes sure that :math:`dH` is zero when
  :math:`Q` is zero. The parameter :math:`\alpha` is automatically chosen such
  that all possible values of :math:`Q` and their accompanying :math:`dH` are
  included.

  With minimization of pumping costs (i.e. power), the optimization results
  will satisfy the first inequality constraint as if it was an equality
  constraint, provided the power relationship of every pump is monotonically
  increasing with `H`.

  .. note::

    Only positive flow is allowed (read: enforced).

  .. cpp:var:: Real C = 0.0


PumpingStation
~~~~~~~~~~~~~~

.. cpp:class:: PumpingStation : Deltares::ChannelFlow::Internal::HQTwoPort

  Represents a pumping station object, containing one or more :cpp:class:`Pump`
  or :cpp:class:`Resistance` objects.

  .. cpp:var:: Integer n_pumps

    The number of pumps contained in the pumping station. This is necessary to
    enforce the right size of e.g. the :cpp:var:`pump_switching_matrix`.

  .. cpp:var:: Integer pump_switching_matrix[n_pumps, n_pumps] = -999

    Together with :cpp:var:`pump_switching_constraints` describes which pumps
    are allowed to be on at the same time. The default value of -999 will make
    Python fill it with the default matrix. This default matrix implies that the
    second pump can only be on when the first pump is on, that the third pump
    can only be on when the second pump is on, etc.

    In matrix multiplication form

    .. math::

      b[:,1] \le A \cdot x \le b[:,2]

    with :math:`A` the :cpp:var:`pump_switching_matrix`, :math:`b` the
    :cpp:var:`pump_switching_constraints`, and :math:`x` the vector of pump statuses:

    .. math::

       x = \left[\begin{array}{cc}S_1 \\ S_2 \\ S_3 \\ \vdots \end{array}\right]

    where :math:`S_1` is the status of pump 1 (on = `1`, off = `0`).

    So the default matrix, where a pump being on requires all lower numbered
    pumps to be on as well, can be expressed as follows:

    .. math::

       A = \left[\begin{array}{ccc}0 & 0 & 0\\1 & -1 & 0\\1 & 1 & -2\end{array}\right]

    with :cpp:var:`pump_switching_constraints` equal to:

    .. math::

       b = \left[\begin{array}{cc}0 & 0\\0 & 1\\0 & 2\end{array}\right]

    To allow all pumps to switch independently from each other, it is sufficient
    to set the coefficient matrix to all zeros (e.g. ``pump_switching_matrix =
    fill(0, n_pumps, n_pumps)``). For rows in the matrix not containing any non-
    zero values, the accompanying constraints are not applied.

    .. note::

      Only square matrices allowed, i.e. a single constraint per pump.

  .. cpp:var:: Integer pump_switching_constraints[n_pumps, 2]

    See discussion in :cpp:var:`pump_switching_matrix`.


Weir
----

.. cpp:namespace:: Deltares::ChannelFlow::Hydraulic::Structures::Weir

.. cpp:class:: Weir : Deltares::ChannelFlow::Hydraulic::Structures

  Represents a general movable-crest weir object described by the conventional
  weir equation (see e.g. Swamee, Prabhata K. "Generalized rectangular weir equations."
  Journal of Hydraulic Engineering 114.8 (1988): 945-949.):

  .. math::

    Q = \frac{2}{3} C B \sqrt{2g} \left( H - H_w \right)^{1.5}

  where `Q` is the discharge of the weir, `C` is the weir discharge
  coefficient (very well approximated by 0.61), `B` is the width of the weir,
  `g` is the acceleration of gravity, `H` is the water level, and :math:`H_w`
  is the level of the movable weir crest. The equation assumes critical flow
  over the weir crest.

  .. cpp:var:: Modelica::SIunits::Length width

    The physical width of the weir.

  .. cpp:var:: Modelica::SIunits::VolumeFlowRate q_min

    The minimal possible discharge on this weir. It can be known from the
    physical characteristics of the system. The linear approximation works the
    best if this is set as tight as possible. It is allowed to set it to zero.

  .. cpp:var:: Modelica::SIunits::VolumeFlowRate q_max

    The maximum physically possible flow over the weir. It should be set as
    tight as possible

  .. cpp:var:: Modelica::SIunits::Length hw_min

    The minimal possible crest elevation.

  .. cpp:var:: Modelica::SIunits::Length hw_max

    The maximum possible crest elevation.

  .. cpp:var:: Real weir_coef = 0.61

    The discharge coefficient of the weir. Typically the default value of
    0.61.

Orifice
-------

.. cpp:namespace:: Deltares::ChannelFlow::Hydraulic::Structures::Orifice

.. cpp:class:: Orifice : Deltares::ChannelFlow::Hydraulic::Structures::Orifice

  
  
  Represents a single (controllable) orifice object. The orifice can only have
  flow in the positive direction, and allows for all flows between zero and a
  head-dependent maximum flow when the head loss is in the right direction.
  This class is wrapped by the Python API's
  :py:class:`~rtctools_hydraulic_structures.orifice_mixin.Orifice` which is
  used to implement the relationship between `H` and `Q` as an inequality
  constraint in Python:

  .. math::

    Q \le C_d \cdot A \cdot \sqrt{2 g \Delta H}

  where :math:`Q` is the pump discharge, and :math:`\Delta H` is the head over
  the orifice (:math:`\Delta H = H_{down} - H_{up}`).

  Note that the convex inequality constraint means that the orifice is
  controllable, i.e. can be opened partially. If the equation above would be
  formulated as an equality constraint, it would result in a non-convex
  problem, which does not fit in the philosophy of RTC-Tools.

  In addition to the inequality constraint above, we also need an integer to
  make sure :math:`Q = 0` in case the the head loss is negative and the
  orifice is closed:

  .. math::

    \begin{aligned} \
    \Delta H - \left(1 - S\right) \cdot M &\le 0\\[5pt]
    \Delta H + S \cdot M &\ge 0\\[5pt]
    0 \le Q + \left(1 - status\right) \cdot Q_{max} &\le Q_{max}\\[5pt]
    \end{aligned}

  where :math:`S` is the status of the orifice (open = 1, closed = 0), and
  :math:`Q_{max}` is the maximum discharge (at the maximum possible head)
  through the orifice.

  .. cpp:var:: Modelica::SIunits::Length dH_max = 10.0

    The maximum possible head over the orifice. Amongst other used to
    calculate the maximum possible discharge through the orifice if it is not
    set using bounds, and to make sure the inequality constraint relating
    :math:`Q` to :math:`H` remains feasible when the orifice is closed.

  .. cpp:var:: Modelica::SIunits::Area area = 1.0

    The cross-sectional area of flow constriction. Note that this is the
    geometric opening area, which is multiplied with the discharge coefficient
    (see below) to get the effective area.

  .. cpp:var:: Real discharge_coefficient = 0.61

    The discharge coefficient, i.e. the ratio of the actual discharge to the
    theoretical discharge.