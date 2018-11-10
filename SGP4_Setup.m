function SGP4_Setup(longstr1, longstr2)
% SGP4_Setup
% Brandon Jackson
% bajackso@mtu.edu
% 9th July 2013
%
%
% Inputs:
% longstr1: TLE character String Line 1
% longstr2: TLE character Sring Line 2
%
% Outputs:
% satrec: structure containing all of the SGP4 satellite information 
%
% Coupling:
% getgravconst
% days2mdhms
% jday
% sgp4init
%
% References:
% Norad Spacetrack Report #3
% Vallado, Crawford, Hujsak, Kelso 2006

%% Define Global Variables
global satrec gravc
% gravc = struct('mu', 0.0,...
%     'radiusearthkm', 0.0,...
%     'xke', 0.0,...
%     'tumin', 0.0,...
%     'j2', 0.0,...
%     'j3', 0.0,...
%     'j4', 0.0,...
%     'j3oj2', 0.0);

%% Include extrinsic functions
% coder.extrinsic('custom_str2double');

%% WGS-72 Earth Constants
% sgp4fix identify constants and allow alternate values
% Options 721 72 84
getgravc( 721 );

%% Define Constants
deg2rad  =   pi / 180.0;         %  0.01745329251994330;  % [deg/rad]
    xpdotp   =  1440.0 / (2.0*pi);   % 229.1831180523293;  % [rev/day]/[rad/min]  

    revnum = 0; 
    elnum  = 0;
    year   = 0; 
    satrec.error = 0;

%% Parse TLE
% Set the implied decimal points since doing a formated read fixes for bad
% input data values (missing, ...)
    for (j = 11:16)
        if (longstr1(j) == ' ')
            longstr1(j) = '_';
        end
    end

    if (longstr1(45) ~= ' ')
        longstr1(44) = longstr1(45);
    end
    longstr1(45) = '.';
     
    if (longstr1(8) == ' ')
        longstr1(8) = 'U';
    end

    if (longstr1(10) == ' ')
        longstr1(10) = '.';
    end

    for (j = 46:50)
        if (longstr1(j) == ' ')
            longstr1(j) = '0';
        end
    end
    if (longstr1(52) == ' ')
        longstr1(52) = '0';
    end
    if (longstr1(54) ~= ' ')
        longstr1(53) = longstr1(54);
    end
    longstr1(54) = '.';

    longstr2(26) = '.';
     
    for (j = 27:33)
        if (longstr2(j) == ' ')
            longstr2(j) = '0';
        end
    end
     
    if (longstr1(63) == ' ')
        longstr1(63) = '0';
    end

    if ((length(longstr1) < 68) || (longstr1(68) == ' '))
        longstr1(68) = '0';
    end

    % parse first line
    carnumb = custom_str2double(longstr1(1));
    satrec.satnum = custom_str2double(longstr1(3:7));
    classification = longstr1(8);
    intldesg = longstr1(10:17);
    satrec.epochyr = custom_str2double(longstr1(19:20));
    satrec.epochdays = custom_str2double(longstr1(21:32));
    satrec.ndot = custom_str2double(longstr1(34:43));
    satrec.nddot = custom_str2double(longstr1(44:50));
    nexp = custom_str2double(longstr1(51:52));
    satrec.bstar = custom_str2double(longstr1(53:59));
    ibexp = custom_str2double(longstr1(60:61));
    numb = custom_str2double(longstr1(63));
    elnum = custom_str2double(longstr1(65:68));
 
    % parse second line
    cardnumb = custom_str2double(longstr2(1));
    satrec.satnum = custom_str2double(longstr2(3:7));
    satrec.inclo = custom_str2double(longstr2(8:16));
    satrec.nodeo = custom_str2double(longstr2(17:25));
    satrec.ecco = custom_str2double(longstr2(26:33));
    satrec.argpo = custom_str2double(longstr2(34:42));
    satrec.mo = custom_str2double(longstr2(43:51));
    satrec.no = custom_str2double(longstr2(52:63));
    revnum = custom_str2double(longstr2(64:68));

% find no, ndot, nddot
    satrec.no   = satrec.no / xpdotp; %//* rad/min
    satrec.nddot= satrec.nddot * 10.0^nexp;
    satrec.bstar= satrec.bstar * 10.0^ibexp;

% convert to sgp4 units
    satrec.a    = (satrec.no*gravc.tumin)^(-2/3);                % [er]
    satrec.ndot = satrec.ndot  / (xpdotp*1440.0);          % [rad/min^2]
    satrec.nddot= satrec.nddot / (xpdotp*1440.0*1440);     % [rad/min^3]

% find standard orbital elements
    satrec.inclo = satrec.inclo  * deg2rad;
    satrec.nodeo = satrec.nodeo * deg2rad;
    satrec.argpo = satrec.argpo  * deg2rad;
    satrec.mo    = satrec.mo     *deg2rad;

    satrec.alta = satrec.a*(1.0 + satrec.ecco) - 1.0;
    satrec.altp = satrec.a*(1.0 - satrec.ecco) - 1.0;

%% Find SGP4 Epoch Time of element set
%  Remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch) and
%  minutes from the epoch (time)

% Temp fix for years 1957-2056
% Correct fix will occur when year is 7-digits in 21e
     if (satrec.epochyr < 57)
         year= satrec.epochyr + 2000;
       else
         year= satrec.epochyr + 1900;
     end;

     [mon,day,hr,minute,sec] = days2mdh ( year,satrec.epochdays );
     satrec.jdsatepoch = jday( year,mon,day,hr,minute,sec );
     
 %% Initialize the orbit at SGP4 Epoch    
 sgp4epoch = satrec.jdsatepoch - 2433281.5; % days since 0 Jan 1950
 sgp4init(sgp4epoch);
end

    function output_double = custom_str2double(input_str)
        i = 0.0;
        output_double = 0;
        npast = 0;
        num_sign = 1;
        input_str = input_str(input_str ~= ' ');
        for i = 1:length(input_str)
            if input_str(i) == '-'
                num_sign = -1; 
            elseif input_str(i) == '.' && npast == 0
                npast = 1;
            elseif npast >= 1;
                output_double = output_double + (input_str(i) - '0')*10^-npast;
                npast = npast+1;
            else
                % Value is greater than 1
                output_double = output_double*10 + (input_str(i) - '0');
            end
        end
        output_double = num_sign*output_double;
    end