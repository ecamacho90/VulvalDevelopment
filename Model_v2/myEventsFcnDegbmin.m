function [value,isterminal,direction]=myEventsFcnDegbmin(t,y,equilibriax,mindistancef1,mindistancef2)

value1 = norm(y-[equilibriax(1);0])-mindistancef1; %It is closer to the point with bigger x

value2 = norm(y-[equilibriax(2);0])-mindistancef2; %It is closer to the point with smaller x . which threshold of closeness should we put?

value = [value1; value2];

isterminal = [1;1];

direction = [0;0];