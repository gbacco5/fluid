function X = SolveSystemMatlabOctave(fun, X0, options, pars)
% SOLVESYSTEMMATLABOCTAVE solves a system through the integrated fsolve
% function.
%
%    SOLVESYSTEMMATLABOCTAVE(fun, X0, options)

% import parameters from main workspace
%names = fieldnames(pars);
%for i=1:length(names)
%  eval([names{i} '=pars.' names{i} ';' ]);
%end


if isOctave
  opt.Display = options.Display; % turn off folve display
  opt.Algorithm = options.Algorithm; % 
  opt.TolFun = options.FunctionTolerance;
  opt.TolX = options.StepTolerance;

  X = fsolve(fun, X0, opt);
  
else
  opt = optimoptions('fsolve');
  opt.Display = options.Display; % turn off folve display
  opt.Algorithm = options.Algorithm; % 
  opt.FunctionTolerance = options.FunctionTolerance;
  opt.StepTolerance = options.StepTolerance;
  
  X = fsolve(fun, X0, opt);
  
  
end