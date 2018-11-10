longstr1 = '1 25544U 98067A   14065.36326021  .00022900  00000-0  39201-3 0  8214';
longstr2 = '2 25544  51.6476 249.2804 0003302 214.0242 215.7852 15.51134586875372';
SGP4_Setup(longstr1, longstr2)

% 4 hours in summer 5 hours in winter do to daylight savings

JD = juliandate(clock + [0 0 0 5 0 0]);
tsince = JD2Tsince(JD);
r = sgp4(tsince);
r = r*1000;
lla = ecef2lla(r');
lla(1:2)

r = [];
for i = 0:0.1:90
    r(end+1,:) = sgp4(i)*1000;
end

figure()
comet3(r(:,1), r(:,2), r(:,3));
