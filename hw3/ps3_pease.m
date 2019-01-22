%% AMATH 481 - PS 3 - Spencer Pease
%

clear all; close all;


%% Exercise 1 - Sparse matrices
%

% Setup -----------------------------------------------------------------------

m = 8;                  % N value in x and y
n = m^2;                % total matrix size
d = (10 - (-10)) / m;   % delta x, delta y

e0 = zeros(n, 1);       % base vector of zeros
e1 = ones(n, 1);        % base vector of ones


% Solve -----------------------------------------------------------------------

%%% Create matrix A %%%

e2 = ones(n, 1); e2(mod(0:n-1, m) == 0) = 0; 
e3 = zeros(n, 1); e3(mod(1:n, m) == 0) = 1;

ind_a = [1 m-1 m (n-m)];
diag_a = [e2 e3 e1 e1];

A = (1 / d^2) * ...
    spdiags([rot90(diag_a, 2) -4*e1 diag_a], [-fliplr(ind_a) 0 ind_a], n, n);


%%% Create matrix B %%%

% ind_b = 0:m:n-1;
% ind_b = ind_b(mod(1:m, 2) == 0);
% 
% tmp_ind = (mod(1:length(ind_b), 2) == 0);
% 
% diag_b = repmat(e1, 1, length(ind_b));
% diag_b(:, tmp_ind) = -1*diag_b(:, tmp_ind);

ind_b = [m n-m];
diag_b = [e1 -e1];

B = (1/(2*d)) * ...
    spdiags([-1*rot90(diag_b, 2), diag_b], [-fliplr(ind_b), ind_b], n, n);


%%% Create matrix C %%%

e4 = repmat([0; ones(m-1, 1)], m, 1);
e5 = repmat([zeros(m-1, 1); -1], m, 1);

ind_c = [1, m-1];
diag_c = [e4, e5];

C = (1/(2*d)) * ...
    spdiags([-1*rot90(diag_c, 2), diag_c], [-fliplr(ind_c), ind_c], n, n);


% Plot and Fit ----------------------------------------------------------------

% imagesc(A)
% imagesc(B)
% imagesc(C)


% Answers ---------------------------------------------------------------------

A1 = full(A);
A2 = full(B);
A3 = full(C);


%% Exercise 2 - FFT 
%

% Setup -----------------------------------------------------------------------

load("Fmat.mat"); load("permvec.mat")


% Solve -----------------------------------------------------------------------

center = Fmat(161:240, 161:240);
chunks = length(permvec);
chunk_size = size(center, 1) / sqrt(chunks);
p_center = zeros(sqrt(chunks) * chunk_size);

[p_row, p_col] = ind2sub([sqrt(chunks), sqrt(chunks)], permvec);
[c_row, c_col] = ind2sub([sqrt(chunks), sqrt(chunks)], 1:chunks);

for j = 1:chunks
    
    cr_ind = (c_row(j)*chunk_size - (chunk_size - 1)):c_row(j)*chunk_size;
    cc_ind = (c_col(j)*chunk_size - (chunk_size - 1)):c_col(j)*chunk_size;
    
    pr_ind = (p_row(j)*chunk_size - (chunk_size - 1)):p_row(j)*chunk_size;
    pc_ind = (p_col(j)*chunk_size - (chunk_size - 1)):p_col(j)*chunk_size;
    
    p_center(cr_ind, cc_ind) = center(pr_ind, pc_ind); 
    
end

decrypt_mat = Fmat; 
decrypt_mat(161:240, 161:240) = p_center;

unshift_mat = ifftshift(decrypt_mat);
transformed = abs(ifft2(unshift_mat));


% Plot and Fit ----------------------------------------------------------------

% set(gcf,'colormap',gray);
% imagesc(uint8(transformed));


% Answers ---------------------------------------------------------------------

A4 = abs(decrypt_mat);
A5 = transformed;


%% Write Data
%

save A1.dat  A1  -ascii
save A2.dat  A2  -ascii
save A3.dat  A3  -ascii
save A4.dat  A4  -ascii
save A5.dat  A5  -ascii

