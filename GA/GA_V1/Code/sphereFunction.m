function [y] = sphereFunction(values)

d = length(values);
sum = 0;

for i = 1:d
	x = values(i);
	sum = sum + x^2;
end

y = sum;

end
