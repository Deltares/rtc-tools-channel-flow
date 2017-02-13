within Deltares.ChannelFlow.Internal;

partial model HQOnePort "Partial model of one port"
  replaceable connector HQPort1 = Deltares.ChannelFlow.Interfaces.HQPort;
  HQPort1 HQ annotation(Placement(visible = true, transformation(origin = {0, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {0, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
end HQOnePort;