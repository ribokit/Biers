function [ensemble, coefficients, optdata] = optimize_RNA_ensemble(sequence, tol, data, use_bonuses)
optdata = zeros(size(data));
expdata = data;
expdata(isnan(data)) = 2;
structnum = 1;
fprintf('Calculating structures ...\n');
if use_bonuses
    structures = get_RNA_structures(sequence, 20, data);
else
    structures = get_RNA_structures(sequence, 20, []);
end
A = zeros(length(sequence), length(structures));
fprintf('Setting constraints\n');
for i=1:length(structures)
    A(:, i) = to_constraints(structures{i});
end
fprintf('Beginning optimization\n');
before = norm(optdata - expdata);
dif = 1;
structures
fprintf('Beginning with ensemble of %d structures \nError %d\n', length(structures), norm(optdata - expdata));
while norm(optdata - expdata) > tol && structnum < length(structures) && dif > 0.001
    fprintf('Number of structures %d\nError %d\n', structnum, norm(optdata - expdata));
    cvx_begin
        variable x(structnum);
        minimize(norm(A(:, 1:structnum)*x - expdata));
        subject to
            x >= 0;
    cvx_end
    optdata = A(:, 1:structnum)*x;
    structnum = structnum + 1;
    dif = abs(norm(optdata - expdata) - before);
    before = norm(optdata - expdata);
end
if exist('x')
    coefficients = x;
    structnum = structnum -1;
else
    coefficients = [1];
end
fprintf('Finished with %d structures\n', structnum);
ensemble = structures(1:structnum);
end

function c = to_constraints(structure)
c = zeros(length(structure), 1);
for i=1:length(structure)
    if structure(i) == '.'
        c(i) = 1;
    end
end
end

