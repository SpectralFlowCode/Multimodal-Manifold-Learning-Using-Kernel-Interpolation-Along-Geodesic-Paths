This is a Matlab code that reproduces the experimental results described in the paper.

To begin the user-interface script execute ``Main.m``.

Execute ``ReprouceAll.m`` in order to reproduce all the results.
In the beginning, the script will download and unzip the datasets.
The generated figure will be saved automatically to a folder named "OutputFigures".
The generated tables will be saved automatically to a folder named "OutputTables".

It takes about ~50 minutes to reproduce all the results.

Execute the following scripts in order to reproduce only part of the results:
- ``Puppets_ReproduceFigures.m``        - to reproduce the figure of the puppets toy example (Figure 1).
- ``Simulations_ReproduceFigures.m``      - to reproduce the figures from the simulations (2DTori - Figures 2+3, 2DFlat - Figures 14+15 in the SM).
- ``ENose_ReproduceFigures.m``            - to reproduce the figures from the  E-Nose experiment (Figure 4 and Figure 9 in the SM).
- ``ENose_ReproduceTable.m``              - to reproduce a table summarizing the objective results on the  E-Nose dataset (Table 1).
- ``CM_ReproduceFigures.m``               - to reproduce the figures from the condition monitoring experiment (Figure 5).
- ``CM_ReproduceTable.m``                 - to reproduce a table summarizing the objective results on the condition monitoring dataset (Table 2).
- ``NoisyMnist_ReproduceFigures.m``       - to reproduce the figures from the noisy mnist experiment (Figure 7).
- ``NoisyMnist_ReproduceTable.m``         - to reproduce a table summarizing the objective results on the noisy mnist dataset (Table 3).
