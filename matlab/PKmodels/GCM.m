function [DSCmodel] = GCM(x,Cp,Ce)
% Gradient Correction Model
% GCM model for DSC. Using Cp and Ce derived using DCE. 

r2vasc = x(1); % r2* of vasculature space
r2ees = x(2);  % r2* of EES space

DSCmodel = (r2vasc*(abs(Cp-Ce)) + r2ees*Ce);
end