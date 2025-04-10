function [r]=fishers_z_to_r(z)

r = (exp(2.*z)-1)./(exp(2.*z)+1);