data = reshape(data, 3, length(data)/3);
i = data(1,:);
j = data(2,:);
s = data(3,:);
H = sparse(i,j,s);
eig_vals = eigs(H);
gnd = min(eig_vals);