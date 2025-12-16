within Deltares.ChannelFlow.HydraulicStructures.PumpingStation;

partial model Pump "Pump with QHP relationship"
  extends Deltares.ChannelFlow.Hydraulic.Structures.Pump;

  // Increasing row number is increasing H power (staring at 0th power).
  // Increasing column number is increasing Q power (staring at 0th power).
  parameter Real power_coefficients[:, :, :];
  parameter Real speed_coefficients[:, :] = {{0.0}};


  // Array of working area polynomials, each a function of Q and H. The
  // coefficients of each polynomial are like the power coefficients, in that
  // increasing row (second index) is increasing power of H, and increasing
  // column (third index) is increasing power of Q.
  parameter Real working_area[:, :, :];

  // For each of the polynomials in the working area we have to specify whether
  // the expression should evaluate to a positive expression (=1), or a
  // negative expression (=-1).
  // NOTE: May become unnecessary to specify this in the future, if we can
  // figure out a way to determine this automatically based on the working
  // area.
  parameter Real working_area_direction[:];

  // Pump's minimum on and off time.
  parameter Modelica.Units.SI.Duration minimum_on = 0.0;
  parameter Modelica.Units.SI.Duration minimum_off = 0.0;

  // NOTE: Enumerations are not supported in JModelica's CasADi interface. We
  // therefore resort to an integer.
  // What head to use for the pump head. This can be
  // -1: The upstream head
  //  0: The differential head (i.e. downstream head minus upstream head)
  //  1: The downstream head.
  parameter Integer head_option = 0;

  parameter Modelica.Units.SI.Energy start_up_energy = 0.0;
  parameter Real start_up_cost = 0.0;

  parameter Modelica.Units.SI.Energy shut_down_energy = 0.0;
  parameter Real shut_down_cost = 0.0;

  // NOTE: The equality constraint setting dH to some combination of HQUp and
  // HQDown (based on head_option) will be added in the Mixin.
  Modelica.Units.SI.Distance dH;

equation
  if head_option == -1 then
    dH = HQUp.H;
  elseif head_option == 1 then
    dH = HQDown.H;
  else
    dH = HQDown.H - HQUp.H;
  end if;
end Pump;
