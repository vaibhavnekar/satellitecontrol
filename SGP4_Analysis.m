% Brandon Jackson
% bajackso@mtu.edu
% 19 August 2013
% This function compares the MATLAB and C coded versions of the SGP4
% propagator against the STK model.

% Two Line Element:
longstr1 = '1 25544U 98067A   13231.43181836  .00009110  00000-0  16523-3 0  4226';
longstr2 = '2 25544  51.6475 154.9600 0004180 335.1118 133.2903 15.50373128844518';

% Start time
StartDate = juliandate([2013 8 19 16 0 0] + [0 0 0 4 0 0]);
PropTimeMin = 2*7*24*60; % two weeks
PropInterval = 1; % min;

TimeArray = 0:PropInterval:PropTimeMin;

global satrec
SGP4_Setup(longstr1, longstr2);

OffSetTimeMin = (StartDate - satrec.jdsatepoch)*24*60;

for indx = 1:length(TimeArray);
    MATLAB_ECEF(1:3,indx) = sgp4(TimeArray(indx) + OffSetTimeMin);
end
MATLAB_ECEF = MATLAB_ECEF';

% Convert STK_LLA = STK_ECEF
[FileName,PathName] = uigetfile('*.csv');
STK_LLA = csvread([PathName FileName]);
STK_LLA = STK_LLA(:,2:4); % truncate rest
STK_LLA(:,3) = STK_LLA(:,3)*1000;
STK_ECEF = lla2ecef(STK_LLA);
STK_ECEF = STK_ECEF/1000;

Error_MATLAB = STK_ECEF - MATLAB_ECEF;
Error_Mag_MATLAB = sqrt(sum(Error_MATLAB.^2,2));

figure();
plot(TimeArray/(24*60), Error_Mag_MATLAB, 'DisplayName', 'MATLAB vs. STK');
legend('-DynamicLegend');
xlabel('Time (days)');
ylabel('Difference (km)');
ylim([0 1]);
grid('on');
