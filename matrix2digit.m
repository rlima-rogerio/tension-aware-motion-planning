%% Encodes the element position in a matrix to a single digit assuming a 
%  continuous counting from left to right, top to bottom.
% INPUT:
%   pos_lin: Row (line) position of the matrix element
%   pos_col: Column position of the matrix element
%   NUM_COL: Number of columns the matrix has
%
% Example:
%   m = [m11 m12 m13        [m1 m2 m3
%        m21 m22 m23    =    m4 m5 m6
%        m31 m32 m33];       m7 m8 m9];
%
% The element m32 is in the linear position 8 
function d = matrix2digit(pos_lin, pos_col, NUM_COL)
    % This edge case accounts for the condition [0,1] = [virtual_start, start]
    if pos_lin == 0
        d = 1;
    else
        d = (pos_lin - 1) * NUM_COL + pos_col;
end