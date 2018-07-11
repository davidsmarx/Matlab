function n = LightMachinerySiIndex(wave)
% n = LightMachinerySiIndex(wave)

% m to um
L = wave/1e-6;

a = 1./(L.^2 - 0.028);
b = 0.013924*a.^2 + 0.138497*a + 3.41696;
c = b - 0.0000209*L.^2;
n = 0.000000148*L.^2 + c;
