%function TwoElementKInterface(epsilon, etaI)
epsilon = 1;
etaI = 1.0;
etaB = 1.0;
K1 = [0 0 0 0; 1 -1 1 -1; -1 1 -1 1; 0 0 0 0];
Keps = K1 - epsilon * transpose(K1);
eta = etaI / 0.5;
K2 = [0 0 0 0; 0 1 -1 0; 0 -1 1 0; 0 0 0 0];
Kinterface = Keps + eta * K2

Kinterior = [2 -2 0 0; -2 2 0 0; 0 0 2 -2; 0 0 -2 2];

eta = etaB / 0.5;
KBl = -[1; 0] * [2 -2] + epsilon * [2; -2] * [1 0] + eta * [ 1 0; 0 0]
KBr = -[0; 1] * [-2 2] + epsilon * [-2; 2] * [0 1] + eta * [ 0 0; 0 1]
Z = zeros(2 ,2);
KB = [KBl Z; Z KBr]

K = Kinterface + Kinterior + KB
a = [0 0.5 0.5 1]';
F = (K * a)'
dK = det(K)
if (abs(dK) < 1e-3)
    [eVec, lambda] = eig(K)
end