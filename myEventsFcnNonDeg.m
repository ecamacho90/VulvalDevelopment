function [value,isterminal,direction]=myEventsFcnNonDeg(t,y,equilibriax,mindistancef1,mindistancef2)

value1 = norm(y-[equilibriax(1);0])-mindistancef1; %It is closer to fate 1

value2 = norm(y-[equilibriax(3);0])-mindistancef2; %It is closer to fate 2

value = [value1; value2];

isterminal = [1;1];

direction = [0;0];