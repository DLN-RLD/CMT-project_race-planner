% GPX to CSV of slopes between points placed every "step" meters 
% Output columns: segment (1 to n-1), slope_percent

% /!\ REMARK : All the local functions are placed at the end of the script for better readability.


%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%

gpxFile = 'track.gpx';
csvFile = 'slopes.csv';
step = 10; % spacing between points [m]

%%%%%%%%%%%%%% 1st step : Read the GPX file %%%%%%%%%%%%%%

[lat, lon, alt] = readGpx(gpxFile);

%%%%%%%%%%%%%% 2nd step : Cumulative horizontal distance %%%%%%%%%%%%%%

[dist, alt] = cumulativeDistance(lat, lon, alt);

%%%%%%%%%%%%%% 3rd step : Resample the altitude every "step" meters %%%%%%%%%%%%%%

s = (0:step:dist(end))'; % s = [0, 10, 20, ...]
altR = interp1(dist, alt, s, 'linear'); % altitude at each of those distances

%%%%%%%%%%%%%% 4th step : Slope between neighbouring points %%%%%%%%%%%%%%

slope = 100 * diff(altR) / step; % slope in percents (n-1 values)
segment = (1:numel(slope))'; %

%%%%%%%%%%%%%% 5th step : Write the CSV %%%%%%%%%%%%%%
Table = table(segment, slope, 'VariableNames', {'segment', 'slope_percent'});
writetable(Table, csvFile);



%%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%

function [lat, lon, alt] = readGpx(gpxFile)
doc = xmlread(gpxFile); % parses the XML file into a tree of Java objects.
nodes = doc.getElementsByTagName('trkpt'); % collects the <trkpt> elements of the XML file.
n = nodes.getLength; % n = number of points taken from the gpx file.
if n < 2
    error('Fewer than 2 <trkpt> found in %s', gpxFile);
end % stops the code in order to prevent errors.
lat = zeros(n,1);
lon = zeros(n,1);
alt = nan(n,1); % NaN, since we don't want a "0" to count as an altitude.
for i = 1:n
    pt = nodes.item(i-1); % "pt" is the current point used in the loop. (the list starts at index 0)
    lat(i) = str2double(char(pt.getAttribute('lat'))); % char() converts the "lat" Java string into a MatLab string. Then str2double converts into a number (not possible directly).
    lon(i) = str2double(char(pt.getAttribute('lon')));
    ele_node = pt.getElementsByTagName('ele'); % can't extract the altitude the same way as "lat" or "lon", because it's a "secondary" attribute (called a "child tag").
    if ele_node.getLength > 0
        alt(i) = str2double(char(ele_node.item(0).getTextContent));
    end
end
alt = fillmissing(alt, 'linear', 'EndValues', 'nearest');  % fills the missing <ele>
end


function [dist, alt] = cumulativeDistance(lat, lon, alt)
R = 6371008.8; % mean Earth radius [m]
phi = deg2rad(lat);
lambda = deg2rad(lon);
dphi = diff(phi); % substracts each latitude to the next one.
dlambda = diff(lambda); % substracts each longitude to the next one.
% The next 2 lines are the "Haversine Formula", used to calculate the
% distance between two points on a sphere.
a = sin(dphi/2).^2 + cos(phi(1:end-1)).*cos(phi(2:end)).*sin(dlambda/2).^2;
seg = 2*R*asin(min(1, sqrt(a))); % segment lengths [m]
dist = [0; cumsum(seg)];
% interp1 needs strictly increasing x: drop zero-length segments
keep = [true; diff(dist) > 0]; % keep is a list of true and false if the condition diff(dist) > 0 is verified
dist = dist(keep); % Only keeps the dist for which the keep value is "true"
alt = alt(keep);
end