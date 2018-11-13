function [Xc] = Crystal()
%% User-defined fractional crystal content
% Crystal content Xc must be between 0 and 1.
Xc = input('What is the crystal fraction? ');

while Xc < 0 || Xc > 1;
    disp('Crystal fraction must be between 0 and 1.');
    % Recall this function if input crystal content is not between 0 and 100%
    [Xc] = Crystal();
end

end

