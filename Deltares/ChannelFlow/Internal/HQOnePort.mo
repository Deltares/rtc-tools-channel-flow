within Deltares.ChannelFlow.Internal;

partial model HQOnePort "Partial model of one port"
  replaceable package medium = Deltares.ChannelFlow.Media.FreshWater;
  Deltares.ChannelFlow.Interfaces.HQCMPort HQ(redeclare package medium = medium) annotation(Placement(visible = true, transformation(origin = {0, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {0, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
end HQOnePort;
