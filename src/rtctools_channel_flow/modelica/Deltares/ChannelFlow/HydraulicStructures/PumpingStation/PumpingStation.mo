within Deltares.ChannelFlow.HydraulicStructures.PumpingStation;

model PumpingStation
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  import SI = Modelica.Units.SI;
  parameter Integer n_pumps = 0;

  // FIXME: For some reason JModelica/CasADi returns {1, 2} for the expression
  // 1:3 if we store it as an Integer, whereas it returns {1, 2, 3} if we
  // store it as a Real. The weird thing is that JModelica does not complain
  // about any size mismatches. Furthermore, transposes also do not seem to
  // work well.
  // To work around these issues, we detect the -999 default array, and
  // overwrite it in Python with the correct one.
  parameter Integer pump_switching_matrix[n_pumps, n_pumps] = fill(-999, n_pumps, n_pumps);
  parameter Integer pump_switching_constraints[n_pumps, 2] = fill(-999, n_pumps, 2);

  SI.VolumeFlowRate Q;
equation
  // Discharge
  Q = HQUp.Q;

  HQUp.M = -HQDown.M;

  // TODO: Annotation / pretty picture. Currently inheriting TwoPort.
end PumpingStation;
