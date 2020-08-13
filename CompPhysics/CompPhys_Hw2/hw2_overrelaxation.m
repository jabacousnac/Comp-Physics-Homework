function hw2_overrelaxation
%%RELAXATION METHOD
format long; %to display to 15 dp
c_vector = [2];
BFN = 100000;%some large number
x_matrix = zeros(BFN, length(c_vector));
x0 = 1;
omega = -1;
for i = 1:length(c_vector)
    for q = 1:BFN
        if q == 1
            x = x0;
        else
            x = (1+omega)*(1 - exp(-c_vector(i)*x)) - omega*x;
            x_matrix(q,i) = x;
            if abs(x_matrix(q,i) - x_matrix(q-1,i)) < 0.000006 %we're done, boys
                n_it = q;
                disp(n_it);
                break;
            end
        end
    end
    out = x_matrix(1:n_it);
    disp(out) %Print out all the values we get from iterations
end
%%

plot((1:length(out)), out, 'r*');
end