function x = gauss_e3(A, B)
    [n, ~] = size(A);
    AB = [A B];  

    for k = 1:n-1
        if AB(k, k) == 0
            for i = k+1:n
                if AB(i, k) ~= 0
                    AB([k i], :) = AB([i k], :);
                    break;
                end
            end
        end
        for i = k+1:n
            if AB(k, k) ~= 0  
                fctr = AB(i, k) / AB(k, k);
                AB(i, k:end) = AB(i, k:end) - fctr * AB(k, k:end);
            end
        end
    end

    x = zeros(n, 1);
    x(n) = AB(n, end) / AB(n, n);
    for i = n-1:-1:1
        x(i) = (AB(i, end) - AB(i, i+1:n) * x(i+1:n)) / AB(i, i);
    end
end

A = [1 2 3; 2 4 7; 3 7 11];
B = [1; 2; 2];
x = gauss_e3(A, B);

disp('Solution:');
disp(x);
disp('The difference Ax-b:');
disp(norm((A*x)-B,2));
