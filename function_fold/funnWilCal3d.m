function [y]=funnWilCal3d(x1,x2,x3)
    y=(1-x1).*cos(pi.*x2)+(x1.^2).*sin(pi.*x3);
end