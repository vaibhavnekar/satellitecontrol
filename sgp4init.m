% -----------------------------------------------------------------------------
%
%                              procedure sgp4init
%
%   this procedure initializes variables for sgp4.
%
% Author:
%   Jeff Beck
%   beckja@alumni.lehigh.edu
%   1.0 (aug 7, 2006) - update for paper dav
% original comments from Vallado C++ version:
%   author        : david vallado                  719-573-2600   28 jun 2005
%
%   inputs        :
%     satn        - satellite number
%     bstar       - sgp4 type drag coefficient              kg/m2er
%     ecco        - eccentricity
%     epoch       - epoch time in days from jan 0, 1950. 0 hr
%     argpo       - argument of perigee (output if ds)
%     inclo       - inclination
%     mo          - mean anomaly (output if ds)
%     no          - mean motion
%     nodeo      - right ascension of ascending node
%
%   outputs       :
%     satrec      - common values for subsequent calls
%     return code - non-zero on error.
%                    1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
%                    2 - mean motion less than 0.0
%                    3 - pert elements, ecc < 0.0  or  ecc > 1.0
%                    4 - semi-latus rectum < 0.0
%                    5 - epoch elements are sub-orbital
%                    6 - satellite has decayed
%
%   locals        :
%     CNODM  , SNODM  , COSIM  , SINIM  , COSOMM , SINOMM
%     Cc1sq  , Cc2    , Cc3
%     Coef   , Coef1
%     cosio4      -
%     day         -
%     dndt        -
%     em          - eccentricity
%     emsq        - eccentricity squared
%     eeta        -
%     etasq       -
%     gam         -
%     argpm       - argument of perigee
%     ndem        -
%     inclm       - inclination
%     mm          - mean anomaly
%     nm          - mean motion
%     perige      - perigee
%     pinvsq      -
%     psisq       -
%     qzms24      -
%     rtemsq      -
%     s1, s2, s3, s4, s5, s6, s7          -
%     sfour       -
%     ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
%     sz1, sz2, sz3
%     sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
%     tc          -
%     temp        -
%     temp1, temp2, temp3       -
%     tsi         -
%     xpidot      -
%     xhdot1      -
%     z1, z2, z3          -
%     z11, z12, z13, z21, z22, z23, z31, z32, z33         -
%
%   coupling      :
%     getgravconst
%     initl       -
%     dscom       -
%     dpper       -
%     dsinit      -
%     sgp4        -
%
%   references    :
%     hoots, roehrich, norad spacetrack report #3 1980
%     hoots, norad spacetrack report #6 1986
%     hoots, schumacher and glover 2004
%     vallado, crawford, hujsak, kelso  2006
%  ----------------------------------------------------------------------------*/

function sgp4init(epoch)
%% Define Global Variables
global satrec gravc

%% Initialization
% Set all near earth variables to zero;

satrec.isimp   = 0;   satrec.method = 'n'; satrec.aycof    = 0.0;
satrec.con41   = 0.0; satrec.cc1    = 0.0; satrec.cc4      = 0.0;
satrec.cc5     = 0.0; satrec.d2     = 0.0; satrec.d3       = 0.0;
satrec.d4      = 0.0; satrec.delmo  = 0.0; satrec.eta      = 0.0;
satrec.argpdot = 0.0; satrec.omgcof = 0.0; satrec.sinmao   = 0.0;
satrec.t       = 0.0; satrec.t2cof  = 0.0; satrec.t3cof    = 0.0;
satrec.t4cof   = 0.0; satrec.t5cof  = 0.0; satrec.x1mth2   = 0.0;
satrec.x7thm1  = 0.0; satrec.mdot   = 0.0; satrec.nodedot = 0.0;
satrec.xlcof   = 0.0; satrec.xmcof  = 0.0; satrec.nodecf  = 0.0;

% Deep space model is not necessary for 700 km orbit which has an expected
% orbiting period of 100ish minutes. These variables were therefore
% deleted.

%% WGS-72 Earth Constants
% sgp4fix identify constants and allow alternate values
% Options 721 72 84
%getgravc( 72 );

ss     = 78.0 / gravc.radiusearthkm + 1.0;
qzms2t = ((120.0 - 78.0) / gravc.radiusearthkm)^4;
x2o3   =  2.0 / 3.0;

% sgp4fix divisor for divide by zero check on inclination
% the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
% 1.5 e-12, so the threshold was changed to 1.5e-12 for consistancy
temp4    =   1.5e-12;

satrec.init = 'y';
satrec.t    = 0.0;

[ainv,  ao, con42, cosio, cosio2, einv, eccsq, omeosq, posq, rp,...
    rteosq, sinio] = initl(epoch);

satrec.error = 0;

% sgp4fix remove this check as it is unnecessary
% the mrt check in sgp4 handles decaying satellite cases even if the starting
%if (rp < 1.0)
%   printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
%    satrec.error = 5;
%end

if ((omeosq >= 0.0 ) || ( satrec.no >= 0.0))
    satrec.isimp = 0;
    if (rp < (220.0 / gravc.radiusearthkm + 1.0))
        satrec.isimp = 1;
    end
    sfour  = ss;
    qzms24 = qzms2t;
    perige = (rp - 1.0) * gravc.radiusearthkm;
    
    % /* - for perigees below 156 km, s and qoms2t are altered - */
    if (perige < 156.0)
        sfour = perige - 78.0;
        if (perige < 98.0)
            sfour = 20.0;
        end
        qzms24 = ((120.0 - sfour) / gravc.radiusearthkm)^4.0;
        sfour  = sfour / gravc.radiusearthkm + 1.0;
    end
    pinvsq = 1.0 / posq;
    
    tsi  = 1.0 / (ao - sfour);
    satrec.eta  = ao * satrec.ecco * tsi;
    etasq = satrec.eta * satrec.eta;
    eeta  = satrec.ecco * satrec.eta;
    psisq = abs(1.0 - etasq);
    coef  = qzms24 * tsi^4.0;
    coef1 = coef / psisq^3.5;
    cc2   = coef1 * satrec.no * (ao * (1.0 + 1.5 * etasq + eeta *...
        (4.0 + etasq)) + 0.375 * gravc.j2 * tsi / psisq * satrec.con41 *...
        (8.0 + 3.0 * etasq * (8.0 + etasq)));
    satrec.cc1   = satrec.bstar * cc2;
    cc3   = 0.0;
    if (satrec.ecco > 1.0e-4)
        cc3 = -2.0 * coef * tsi * gravc.j3oj2 * satrec.no * sinio / satrec.ecco;
    end
    satrec.x1mth2 = 1.0 - cosio2;
    satrec.cc4    = 2.0* satrec.no * coef1 * ao * omeosq *...
        (satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecco *...
        (0.5 + 2.0 * etasq) - gravc.j2 * tsi / (ao * psisq) *...
        (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq *...
        (1.5 - 0.5 * eeta)) + 0.75 * satrec.x1mth2 *...
        (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * satrec.argpo)));
    satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *...
        (etasq + eeta) + eeta * etasq);
    cosio4 = cosio2 * cosio2;
    temp1  = 1.5 * gravc.j2 * pinvsq * satrec.no;
    temp2  = 0.5 * temp1 * gravc.j2 * pinvsq;
    temp3  = -0.46875 * gravc.j4 * pinvsq * pinvsq * satrec.no;
    satrec.mdot     = satrec.no + 0.5 * temp1 * rteosq * satrec.con41 +...
        0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
    satrec.argpdot  = -0.5 * temp1 * con42 + 0.0625 * temp2 *...
        (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +...
        temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
    xhdot1            = -temp1 * cosio;
    satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +...
        2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
    xpidot            =  satrec.argpdot+ satrec.nodedot;
    satrec.omgcof   = satrec.bstar * cc3 * cos(satrec.argpo);
    satrec.xmcof    = 0.0;
    if (satrec.ecco > 1.0e-4)
        satrec.xmcof = -x2o3 * coef * satrec.bstar / eeta;
    end
    satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1;
    satrec.t2cof   = 1.5 * satrec.cc1;
    
    % // sgp4fix for divide by zero with xinco = 180 deg
    if (abs(cosio+1.0) > 1.5e-12)
        satrec.xlcof   = -0.25 * gravc.j3oj2 * sinio *...
            (3.0 + 5.0 * cosio) / (1.0 + cosio);
    else
        satrec.xlcof   = -0.25 * gravc.j3oj2 * sinio *...
            (3.0 + 5.0 * cosio) / temp4;
    end
    satrec.aycof   = -0.5 * gravc.j3oj2 * sinio;
    satrec.delmo   = (1.0 + satrec.eta * cos(satrec.mo))^3;
    satrec.sinmao  = sin(satrec.mo);
    satrec.x7thm1  = 7.0 * cosio2 - 1.0;
    
 % Deep space inititiation removed since satellite orbit period will not
 % require a deep space model.
 
 % Set variables for non deep space model.
    if (satrec.isimp ~= 1)
        cc1sq = satrec.cc1 * satrec.cc1;
        satrec.d2 = 4.0 * ao * tsi * cc1sq;
        temp = satrec.d2 * tsi * satrec.cc1 / 3.0;
        satrec.d3 = (17.0 * ao + sfour) * temp;
        satrec.d4 = 0.5 * temp * ao * tsi *...
            (221.0 * ao + 31.0 * sfour) * satrec.cc1;
        satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
        satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 *...
            (12.0 * satrec.d2 + 10.0 * cc1sq));
        satrec.t5cof = 0.2 * (3.0 * satrec.d4 +...
            12.0 * satrec.cc1 * satrec.d3 +...
            6.0 * satrec.d2 * satrec.d2 +...
            15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));
    end
end % // if omeosq = 0 ...

% Finally propogate to zero epoch to initialise all others.
% sgp4fix take out check to let satellites process until they are actually
% below earth surface
%if(satrec.error == 0)

t_0 = 0.0;
[~] = sgp4(t_0);
%end

satrec.init = 'n';

% global idebug dbgfile
% if idebug
%     debug6;
% end

return;
