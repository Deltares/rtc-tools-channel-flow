within Deltares.ChannelFlow.Internal;

partial model HQOnePort "Partial model of one port"
  replaceable connector HQPort = Deltares.ChannelFlow.Interfaces.HQPort;
  HQPort HQ annotation(Placement(visible = true, transformation(origin = {0, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {0, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
end HQOnePort;