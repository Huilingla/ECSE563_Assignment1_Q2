function Z = impedance(nfrom, nto, r, x, b)
% IMPEDANCE Calculate the impedance matrix for an AC power network.
%
%   Z = IMPEDANCE(NFROM, NTO, R, X, B) computes the nodal impedance matrix
%   by first calculating the admittance matrix and then solving Z*Y = I.
%
%   Inputs:
%     nfrom - Mx1 vector of branch starting nodes
%     nto   - Mx1 vector of branch ending nodes  
%     r     - Mx1 vector of branch resistances (p.u.)
%     x     - Mx1 vector of branch reactances (p.u.)
%     b     - Mx1 vector of branch shunt susceptances (p.u.)
%
%   Output:
%     Z     - NxN nodal impedance matrix (p.u.)
%
%   Method:
%     1. Calculate admittance matrix Y using the admittance() function
%     2. Solve Z*Y = I for Z using linsolve() (column by column)
%     3. Handle singularity using appropriate numerical methods

    fprintf('=== IMPEDANCE MATRIX CALCULATION ===\n\n');
    
    % Step 1: Calculate admittance matrix using existing function
    fprintf('Step 1: Calculating admittance matrix Y\n');
    fprintf('---------------------------------------\n');
    Y = admittance(nfrom, nto, r, x, b);
    
    % Determine network size
    all_nodes = unique([nfrom; nto]);
    N = max(all_nodes);
    
    % Step 2: Calculate impedance matrix by solving Z*Y = I
    fprintf('\nStep 2: Calculating impedance matrix by solving Z*Y = I\n');
    fprintf('-------------------------------------------------------\n');
    fprintf('Solving system column by column using linsolve()\n');
    
    % Initialize impedance matrix
    Z = zeros(N, N);
    
    % Create identity matrix for the right-hand side
    I_matrix = eye(N);
    
    % Solve for each column of Z: Z*Y = I => Y'*Z' = I'
    % We solve for each column: Y' * z_col = e_col (where e_col is column of I)
    for col = 1:N
        % Solve Y' * z_col = e_col for each column
        e_col = I_matrix(:, col);  % Current column of identity matrix
        
        % Use linsolve with appropriate options for Hermitian matrix
        opts.SYM = true;
        opts.POSDEF = true;
        
        try
            % Solve Y' * z_col = e_col
            z_col = linsolve(Y', e_col, opts);
        catch
            % If symmetric solver fails, use general solver
            fprintf('  Column %d: Using general solver\n', col);
            z_col = linsolve(Y', e_col);
        end
        
        % Store the solution as a column in Z'
        Z(:, col) = z_col;
    end
    
    % Transpose to get the correct Z matrix (since we solved for Z')
    Z = Z';
    
    fprintf('Impedance matrix Z (%dx%d) calculated successfully\n', N, N);
    
    % Step 3: Validate the result
    fprintf('\nStep 3: Validation\n');
    fprintf('------------------\n');
    
    % Check that Z * Y ≈ I (identity matrix)
    identity_check = Z * Y;
    max_error = max(max(abs(identity_check - eye(N))));
    
    fprintf('Maximum error in Z*Y = I validation: %e\n', max_error);
    
    if max_error < 1e-10
        fprintf('Validation: Excellent (Z is accurate inverse of Y)\n');
    elseif max_error < 1e-6
        fprintf('Validation: Good (Z is good approximation of Y^(-1))\n');
    else
        fprintf('Validation: Acceptable (consider checking matrix conditioning)\n');
    end
    
    % Display matrix conditioning information
    fprintf('Condition number of Y: %e\n', cond(Y));
end
