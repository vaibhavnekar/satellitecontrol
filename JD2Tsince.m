function tsince = JD2Tsince(JD)
% JD2Tsince
% Brandon Jackson
% bajackso@mtu.edu
% 9th July 2013
%
% This function uses the julian date of the epoch stored in the satrec
% structure to determine the time since the epoch given the current julian
% date.

global satrec

tsince = (JD - satrec.jdsatepoch)*24*60;
