# FRL_Income_Inequality

Francisco Tavares Garcia and Jamie L. Cross (2024). 
The impact of monetary policy on income inequality: Does inflation targeting matter?
Finance Research Letters, 2024, 61, 105006.
https://doi.org/10.1016/j.frl.2024.105006.

-----------------------------------------------------------------------------------------------------------------------------------------------
REPLICATION FILES
-----------------------------------------------------------------------------------------------------------------------------------------------
The replication folder has 3 subfolders:
1. Data: contains fullset.xlsx, which is the full sample for the G7 economies. For each country, there are 3 subsamples separated by sheet.
2. Functions: contains various sub-functions used in the main run file.
3. Results: A storage folder where the main run file will save the Impulse Response Functions (IRF) MATLAB files.

The replication codes provided are in MATLAB. The main script is "Fig3_by_country.m". Running this script will load the data from “fullset.xlsx” from the data folder, select the country of choice, and generate the impulse responses for each of the country's subsamples. The IRF figure files are saved in the folder “results/irf”.
