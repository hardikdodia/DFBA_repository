function combinations = MakeCombinations(n)
    x = floor(0.8*n);
    combinations = zeros(n,1);
    combinations(1:x) = 2;
    combinations(x+1:end) = 1;
    combinations = combinations.';
end