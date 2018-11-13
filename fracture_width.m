function [mean_wf] = fracture_width(nf,l)
% Ensures that cumulative fracture width is not greater than the with of
% intact host material.

mean_wf = input('What is the mean fracture width [m]? ');
while mean_wf*nf >= l*0.5;
    disp('Cumulative fracture width cannot be greater than intact material.');
    [mean_wf] = fracture_width(nf,l);
end
end

