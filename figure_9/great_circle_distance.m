function d = great_circle_distance(az1,el1,az2,el2,r)
% INPUTS
% az1, az2, el1, ez2 : azimuth and elevation coordinate vectors of particle 1 and 2, in radians
% r      : radius of sphere

% Using vincenty formula for accuracy:
% 1975, Vincenty T. "Direct and Inverse Solutions of Geodesics
%         on the Ellipsoid with Application of Nested Equations."
%         Survey Review. Directorate of Overseas Surveys. 23(176):88-93
% Note that the azimuth (az) is the same as Vincenty's longitude (lambda)
%             elevation (el) is the same as Vincenty's latitude (phi)
sigma = atan( sqrt( (cos(el2).*sin(abs(az1-az2))).^2 +...
                    (cos(el1).*sin(el2)-sin(el1).*cos(el2).*cos(abs(az1-az2))).^2)  / ...
            (sin(el1).*sin(el2)+cos(el1).*cos(el2).*cos(abs(az1-az2))));
d=sigma*r;
end