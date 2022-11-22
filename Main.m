disp('--------------------------------------------------------------------------------------------------------');
disp('This is a Matlab code that reproduces the experimental results described in the paper.')
disp('In the beginning, the script will download and unzip the datasets.')
disp('The generated figure will be saved automatically to a folder named "./OutputFigures/".')
disp('The generated tables will be saved automatically to a folder named "./OutputTables/".')

disp('Type:')
disp(' - "All"                   - to reproduce all the results (takes about ~50 minutes)')
disp(' - "Puppets"               - to reproduce the figure of the puppets toy example (Figure 1)')
disp(' - "Simulations"           - to reproduce the figures from the simulations (2DTori - Figures 2+3, 2DFlat - Figures 14+15 in the SM).')
disp(' - "ENose-Figures"         - to reproduce the figures from the  E-Nose experiment (Figure 4 and Figure 9 in the SM)')
disp(' - "ENose-Table"           - to reproduce a table summarizing the objective results on the  E-Nose dataset (Table 1).')
disp(' - "CM-Figures"            - to reproduce the figures from the condition monitoring experiment (Figure 5)..')
disp(' - "CM-Table"              - to reproduce a table summarizing the objective results on the condition monitoring dataset (Table 2).')
disp(' - "NoisyMnist-Figures"    - to reproduce the figures from the noisy mnist experiment (Figure 7)')
disp(' - "NoisyMnist-Table"      - to reproduce a table summarizing the objective results on the noisy mnist dataset (Table 3).')
prompt = "Please type your choice:\n";
input_str = input(prompt,'s');
while ~any(strcmp({'all','puppets','simulations','enose-figures','cm-figures','cm-table','noisymnist-figures','noisymnist-table'},lower(input_str)))
    input_str = input( "Not a valid input, please type again:\n",'s');
end
    
switch  lower(input_str)
    case lower('All')
        ReproduceAll()
    case lower('Puppets')
        Puppets_ReproduceFigures()
    case lower('Simulations')
        Simulations_ReproduceFigures()
    case lower('ENose-Figures')
        ENose_ReproduceFigures()
    case lower('ENose-Table')
        ENose_ReproduceTable()
    case lower('CM-Figures')
        CM_ReproduceFigures()
    case lower('CM-Table')
        CM_ReproduceTable
    case lower('NoisyMnist-Figures')
        NoisyMnist_ReproduceFigures()
    case lower('NoisyMnist-Table')
        NoisyMnist_ReproduceTable()
end

