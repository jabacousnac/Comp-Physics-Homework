function hw2_wien
%%
% Function is 5e^-x + x - 5. We see 0 is a root. But we want to find the
% other root
xl = 0.1;
xr = 10;
int = 1;
target = 0.000001;
while int
    left = 5*exp(-xl)+xl-5;
    %right = 5*exp(-xr)+xr-5;
    x = (xl+xr)/2;
    center = 5*exp(-x)+x-5;
    if abs(xr - xl)<target %we're donzo
        int = 0;
        disp(x);
    else
        if left*center > 0;
            xl = x;
        else
            xr = x;
        end
    end
end