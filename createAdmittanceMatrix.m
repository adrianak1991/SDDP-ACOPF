%%% Function returns the admittance matrix and its transformation required
%%% to solve the OPF problem based on the network information.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[admittance] = createAdmittanceMatrix(data_struct)

n = data_struct.net.n;
l = data_struct.net.l;
L = data_struct.net.L;
base = data_struct.net.base;

branch_info =  data_struct.sets.branch_info;
bus_info = data_struct.sets.bus_info;

%% CREATING AN ADMITTANCE MATRIX
Y = zeros(n,n);
y_ff = zeros(l,1);
y_ft = zeros(l,1);
y_tf = zeros(l,1);
y_tt = zeros(l,1);

% definition of additional parameters
status = branch_info(:,11);

ratio = ones(l, 1);                              %% default tap ratio = 1
i = find(branch_info(:, 9));                       %% indices of non-zero tap ratios
ratio(i) = branch_info(i, 9);
shift = branch_info(:,10);

tap= ratio .* exp(1j*pi/180 * shift);

y_sa = status./(branch_info(:,3)+1j*branch_info(:,4));
y_b = status.*branch_info(:,5);

y_tt = y_sa+1j*y_b/2;
y_ff = y_tt./(tap.*conj(tap));
y_ft = -y_sa./conj(tap);
y_tf = -y_sa./tap;

% shunt admittances associated with buses
Y_sh = (bus_info(:,5)+1j*bus_info(:,6))/base;

f = branch_info(:,1);
t = branch_info(:,2);

Y = sparse([f;f;t;t], [f;t;f;t], [y_ff;y_ft;y_tf;y_tt], n, n) + sparse(1:n, 1:n, Y_sh, n, n);
G = real(Y);
B = imag(Y);

admittance.Y = Y;
admittance.G = G;
admittance.B = B;

%% CREATING MATRICES USED IN REFORMULATION
% Last coordinate goes through all buses/lines.
bY_k = zeros(n,2*n,2*n); % for every bus we have different matrix
hbY_k = zeros(n,2*n,2*n); % for every bus we have different matrix
M_k = zeros(n,2*n,2*n); % for every bus we have different matrix

bY_lm = zeros(2*l,2*n,2*n); % for every line we have different matrix (we inlcude both directions)
hbY_lm = zeros(2*l,2*n,2*n); % for every line we have different matrix (we inlcude both directions)
M_lm = zeros(l,2*n,2*n); % for every line we have different matrix (we inlcude both directions)

for i=1:n
    e = zeros(n,1);
    e(i) = 1;
    
    Y_k = e*e'*Y;
    bY_k(i,:,:) = [real(Y_k+transpose(Y_k)) imag(-Y_k+transpose(Y_k));...
        imag(Y_k-transpose(Y_k)) real(Y_k+transpose(Y_k))]/2;
    hbY_k(i,:,:) = -[imag(Y_k+transpose(Y_k)) real(Y_k-transpose(Y_k));...
        real(-Y_k+transpose(Y_k)) imag(Y_k+transpose(Y_k))]/2;
    
    M_k(i,:,:) = [e*e' zeros(n,n); zeros(n,n) e*e'];
end

for i=1:l
    e_l = zeros(n,1);
    e_m = zeros(n,1);
    e_l(L(i,1)) = 1;
    e_m(L(i,2)) = 1;
    
    Y_lm = y_ff(i)*e_l*e_l'+y_ft(i)*e_l*e_m'; 
    bY_lm(i,:,:) = [real(Y_lm+transpose(Y_lm)) imag(-Y_lm+transpose(Y_lm));...
        imag(Y_lm-transpose(Y_lm)) real(Y_lm+transpose(Y_lm))]/2;
    hbY_lm(i,:,:) = -[imag(Y_lm+transpose(Y_lm)) real(Y_lm-transpose(Y_lm));...
        real(-Y_lm+transpose(Y_lm)) imag(Y_lm+transpose(Y_lm))]/2;
    
    Y_lm = y_tt(i)*e_m*e_m'+y_tf(i)*e_m*e_l'; 
    bY_lm(i+l,:,:) = [real(Y_lm+transpose(Y_lm)) imag(-Y_lm+transpose(Y_lm));...
        imag(Y_lm-transpose(Y_lm)) real(Y_lm+transpose(Y_lm))]/2;
    hbY_lm(i+l,:,:) = -[imag(Y_lm+transpose(Y_lm)) real(Y_lm-transpose(Y_lm));...
        real(-Y_lm+transpose(Y_lm)) imag(Y_lm+transpose(Y_lm))]/2;
    
    M_lm(i,:,:) = [(e_l-e_m)*(e_l-e_m)' zeros(n,n); zeros(n,n) (e_l-e_m)*(e_l-e_m)'];
end

admittance.bY_k = bY_k; 
admittance.hbY_k = hbY_k;
admittance.M_k = M_k;

admittance.bY_lm = bY_lm;
admittance.hbY_lm = hbY_lm;
admittance.M_lm = M_lm;

end