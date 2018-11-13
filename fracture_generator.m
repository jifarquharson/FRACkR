function [w] = fracture_generator(window_min, window_max,nf, mean_wf)
%% Generates fracture depths up to a given density
%  Fractures have a equal likelihood of being generated at any point
%  within the user-specified depth range. If the distance between two
%  fractures is less than or equal to the mean fracture width, the
%  function is re-run until this is not the case.
pd = makedist('Uniform', 'lower',window_min,'upper',window_max);
w = random(pd,nf,1); w = sort(w(:));
diffs = diff(w);
if any(diffs)<= (mean_wf/2);
    [w] = fracture_generator(window_min, window_max,nf, mean_wf);
end
end


