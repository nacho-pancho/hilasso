% 
% draw samples from a Jeffreys of Exponential (JOE) distribution
% 
% inputs: 
% N ........ number of samples to draw
% a  ....... start of range for Jeffreys
% b  ....... end of range.
%
% outputs:
% x ........ a vector of N samples from this distribution.
%
% author: Ignacio Ramírez ($Author$)
% date  : $Date$
% rev   : $Revision$
%
function x=mdlsJoeSample(N,a,b)
%
% use the X=F^-1(A) technique with A chosen uniformly between 0
% and 1. Since F cannot be computed in closed form, we need a Lookup table
% and then we interpolate linearly
%
% The lookup table goes from 0 to b*10
%
