function [ y, grad ] = func( x )
%FUNC Summary of this function goes here
%   Detailed explanation goes here
y = (x - 3)*(x - 3);

grad = 2*(x - 3);

end

