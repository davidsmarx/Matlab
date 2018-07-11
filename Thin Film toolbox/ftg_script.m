% Active X example using Film Star Design

% First open Server
Design = actxserver('Ftgdesign1.clsBasic');

% list valid properties:
Design.Design = '0.25H 0.25L 0.25H';
invoke(Design,'CalcPlot');

lam = Design.Spectrum_X;
spectrum = num2cell(double(Design.Spectrum_Y),1);
[Rphi, R, Tphi, T] = deal(spectrum{:});

figure, plot(lam,[Rphi Tphi]/pi), grid
figure, plot(lam,[R T]), grid, legend('R','T')

value = invoke(Design,'IndexCoeff',3,1)

release(Design);
