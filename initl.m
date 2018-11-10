% -----------------------------------------------------------------------------
%
%                            procedure initl
%
%   this procedure initializes the spg4 propagator. all the initialization is
%     consolidated here instead of having multiple loops inside other routines.
%
% Author: 
%   Jeff Beck 
%   beckja@alumni.lehigh.edu
%   1.0 (aug 7, 2006) - update for paper dav
%   1.1 nov 16, 2007 - update for better compliance
% original comments from Vallado C++ version:
%   author        : david vallado                  719-573-2600   28 jun 2005
%
%   inputs        :
%     satrec.ecco        - eccentricity                           0.0 - 1.0
%     epoch       - epoch time in days from jan 0, 1950. 0 hr
%     satrec.inclo       - inclination of satellite
%     satrec.no          - mean motion of satellite
%     satrec.satnum        - satellite number
%
%   outputs       :
%     ainv        - 1.0 / a
%     ao          - semi major axis
%     satrec.con41       -
%     con42       - 1.0 - 5.0 cos(i)
%     cosio       - cosine of inclination
%     cosio2      - cosio squared
%     einv        - 1.0 / e
%     eccsq       - eccentricity squared
%     satrec.method      - flag for deep space                    'd', 'n'
%     omeosq      - 1.0 - satrec.ecco * satrec.ecco
%     posq        - semi-parameter squared
%     rp          - radius of perigee
%     rteosq      - square root of (1.0 - satrec.ecco*satrec.ecco)
%     sinio       - sine of inclination
%     satrec.gsto        - gst at time of observation               rad
%     satrec.no          - mean motion of satellite
%
%   locals        :
%     ak          -
%     d1          -
%     del         -
%     adel        -
%     po          -
%
%   coupling      :
%     gstime      - find greenwich sidereal time from the julian date
%
%   references    :
%     hoots, roehrich, norad spacetrack report #3 1980
%     hoots, norad spacetrack report #6 1986
%     hoots, schumacher and glover 2004
%     vallado, crawford, hujsak, kelso  2006
%  ----------------------------------------------------------------------------*/

function [ainv,  ao, con42, cosio, cosio2, einv, eccsq, omeosq, posq, rp,...
    rteosq, sinio] = initl(epoch);

   % /* -------------------- wgs-72 earth constants ----------------- */
   %     // sgp4fix identify constants and allow alternate values
   x2o3   = 2.0 / 3.0;
   global opsmode
   global satrec gravc

   % /* ------------- calculate auxillary epoch quantities ---------- */
   eccsq  = satrec.ecco * satrec.ecco;
   omeosq = 1.0 - eccsq;
   rteosq = sqrt(omeosq);
   cosio  = cos(satrec.inclo);
   cosio2 = cosio * cosio;

   % /* ------------------ un-kozai the mean motion ----------------- */
   ak    = (gravc.xke / satrec.no)^x2o3;
   d1    = 0.75 * gravc.j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
   del   = d1 / (ak * ak);
   adel  = ak * (1.0 - del * del - del *...
       (1.0 / 3.0 + 134.0 * del * del / 81.0));
   del   = d1/(adel * adel);
   satrec.no    = satrec.no / (1.0 + del);

   ao    = (gravc.xke / satrec.no)^x2o3;
   sinio = sin(satrec.inclo);
   po    = ao * omeosq;
   con42 = 1.0 - 5.0 * cosio2;
   satrec.con41 = -con42-cosio2-cosio2;
   ainv  = 1.0 / ao;
   einv  = 1.0 / satrec.ecco;
   posq  = po * po;
   rp    = ao * (1.0 - satrec.ecco);
   satrec.method = 'n';

   % sgp4fix modern approach to finding sidereal time
   satrec.gsto = gstime(satrec.jdsatepoch);

%    global idebug dbgfile
%    if isempty(idebug)
%        idebug = 0;
%    elseif idebug
%        debug5;
%    end

   return;
