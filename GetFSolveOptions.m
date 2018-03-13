function opt = GetFSolveOptions(opt)

if isOctave
  opt.Display = opt.Display; % turn off folve display
  opt.Algorithm = opt.Algorithm; %
  opt.TolFun = opt.FunctionTolerance;
  opt.TolX = opt.StepTolerance;
  
end

end
