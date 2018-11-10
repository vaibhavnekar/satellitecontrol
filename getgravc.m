% -----------------------------------------------------------------------------
%
%                           function getgravc
%
%  this function gets constants for the propagator. note that gravc.mu is identified to
%    facilitiate comparisons with newer models.
%
%  author        : david vallado                  719-573-2600   21 jul 2006
%
%  inputs        :
%    whichconst  - which set of constants to use  721, 72, 84
%
%  outputs       :
%    tumin       - minutes in one time unit
%    mu          - earth gravitational parameter
%    radiusearthkm - radius of the earth in km
%    xke         - reciprocal of gravc.tumin
%    j2, j3, j4  - un-normalized zonal harmonic values
%    j3oj2       - j3 divided by j2
%
%  locals        :
%
%  coupling      :
%
%  references    :
%    norad spacetrack report #3
%    vallado, crawford, hujsak, kelso  2006
% [gravc.tumin, gravc.mu, gravc.radiusearthkm, gravc.xke, j2, j3, j4, j3oj2] = getgravc(whichconst);
%  --------------------------------------------------------------------------- */

function getgravc(whichconst)

	global gravc
    
      switch whichconst
          case 721
           % -- wgs-72 low precision str#3 constants --
           gravc.mu     = 398600.79964;        %// in km3 / s2
           gravc.radiusearthkm = 6378.135;     %// km
           gravc.xke    = 0.0743669161;
           gravc.tumin  = 1.0 / gravc.xke;
           gravc.j2     =   0.001082616;
           gravc.j3     =  -0.00000253881;
           gravc.j4     =  -0.00000165597;
           gravc.j3oj2  =  gravc.j3 / gravc.j2;
          case 72
           % ------------ wgs-72 constants ------------
           gravc.mu     = 398600.8;            %// in km3 / s2
           gravc.radiusearthkm = 6378.135;     %// km
           gravc.xke    = 60.0 / sqrt(gravc.radiusearthkm*gravc.radiusearthkm*gravc.radiusearthkm/gravc.mu);
           gravc.tumin  = 1.0 / gravc.xke;
           gravc.j2     =   0.001082616;
           gravc.j3     =  -0.00000253881;
           gravc.j4     =  -0.00000165597;
           gravc.j3oj2  =  gravc.j3 / gravc.j2;
          case 84
           % ------------ wgs-84 constants ------------
           gravc.mu     = 398600.5;            %// in km3 / s2
           gravc.radiusearthkm = 6378.137;     %// km
           gravc.xke    = 60.0 / sqrt(gravc.radiusearthkm*gravc.radiusearthkm*gravc.radiusearthkm/gravc.mu);
           gravc.tumin  = 1.0 / gravc.xke;
           gravc.j2     =   0.00108262998905;
           gravc.j3     =  -0.00000253215306;
           gravc.j4     =  -0.00000161098761;
           gravc.j3oj2  =  gravc.j3 / gravc.j2;
      end;  % case

