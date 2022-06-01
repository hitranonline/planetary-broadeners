# planetary-broadeners
==============================================================================
License

See the LICENSE.txt file available in the "Broadening .py Files" folder.

These HITRAN Broadening Python files are distributed under the Apache 2.0 license

Copyright 2022 HITRAN Team. You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0

==============================================================================

Citation

Y. Tan, F.M. Skinner, S. Samuels, R.J. Hargreaves, R. Hashemi, I.E. Gordon (2022), Submitted to the Astrophysical Journal Supplement Series in March 2022
"H2, He, and CO2 pressure-induced parameters for the HITRAN database. Part II: Line lists of CO2, N2O, CO, SO2, OH, OCS, H2CO, HCN, PH3, H2S and GeH4"

==============================================================================

Purpose

The broadening Python files are designed to apply broadening parameters to HITRAN line lists through dependency on rotational quantum numbers (J+0.2Ka or |m| values).

==============================================================================

Explanation of the Input & Output Files

Input files are given in the "Input Broadening Files" folder and are sampled portions of the HITRAN2020 line list data (not the complete line list).
In particular, the input files are HITRAN line lists in default format called .par (see https://hitran.org/docs/definitions-and-units/ for an explanation of the .par format)

We highly recommend that users directly download the current line lists from the HITRAN database (https://hitran.org/) prior to applying these broadening codes.

These broadening codes can be used on any line list, however the user must remember to alter the format of the read-in file if they are not using a HITRAN .par formatted line list.
If the user is not using a HITRAN .par line list, then the primary columns to isolate when formatting the read-in line list are the rotational quantum number(s) and the branch columns 
(for some molecules the branches are not needed). This is due to the fact that the main data columns used to apply broadening parameters are dependent on rotational quantum number
(J lower, and for some molecules Ka lower) and Branch letter (P, Q or R).

Output files are available in the "Output Broadening Files" folder. The output files are examples of what users should receive if they run the given input files correctly 
from the "Input Broadening Files" folder by using the broadening Python scripts, given in the "Broadening .py Files" folder.

==============================================================================

Accessing Broadening Parameters via HITRANonline

To access these additional foreign-broadening parameters via the HITRAN database directly, 
HITRAN users can proceed to https://hitran.org/ to create a customized output file when downloading line-by-line data.

Here is an example of adding H2 broadening parameters to a line list output within the HITRAN database.
The steps are given as follows:
	Login (only requires your name and email address)
	Go to Line-by-Line under Data Access
	Choose a molecule
	Choose the isotopologues of the molecule
	Specify your spectral range
	(!STOP! This step is the important part where HITRAN is asking what the output format should be)
	Scroll down and click the green button (left-hand side) that says "Create New Output Format"
	(Now you should see a page that says "New Output Format")
	Under "Available Parameters" choose ".par line"
	Then under "Available Parameters" choose "γH2" then "ηH2" and "δH2"
	Give your output format a name within the "Output Format Name:" box
	Give your output format a description within the "Description:" box
	(You can also change the field separator and line endings, also you can apply a fixed width format and an output header line)
	Now click the green button "Save and Return to Data Search"
	(You might get a warning that reads:
	"You have not selected "Global isotopologue ID" or "Molecule ID" and "Isotopologue ID" 
	among of your parameters and so may have no way to know which transition corresponds to which species. Proceed?" 
	Do not worry about this warning, we have the both the "Molecule ID" and "Isotopologue ID" within the HITRAN default .par format, 
	which we already added to our output. Therefore, please proceed with confidence even if you receive this warning.)
	Press "Okay" to exit the warning
	Press "Start Data Search" to retrieve the final line list and its corresponding references with the broadening parameters included.

==============================================================================

Definitions

|m| values are the rotational running index. |m| is defined with the following relationships to the J rotational quanta:
						P-branch: m = -J"
						Q-branch: m = J"
						R-branch: m = J"+1

Broadening parameters refer to the following:
	1) Lorentzian half-widths at half maximum (HWHM) denoted as γH2, γHe, γCO2 and γH2O for H2-, He-, CO2- and H2O-broadening respectively.
	2) Temperature dependence (exponent) of these half widths, denoted as ηH2, ηHe, ηCO2, ηH2O and defined through the power law (see below).
	3) Collisional line shifts, denoted as δH2, δHe, δCO2, δH2O which at the moment are available only for some HITRAN molecules 
	   and in some cases only for some specifically selected lines of these molecules. 

The power law equation for determining the HWHM at T (Temperature) is given as: γ(T) = γ(T0)([T0/T]^η)
	where T0 is the reference temperature (296K in HITRAN) and γ(T0) is the HWHM at the reference temperature.

The broadening Python files utilize 3rd-to-4th order Padé approximants (equation given below) to populate the broadening data throughout the line list.

The Padé approximants were fit to available laboratory data sets in order to extrapolate to other transitions where data was not available.

The 3rd-to-4th order Padé approximant: γx(|m|) = (a0 + a1|m| + a2|m|^2 + a3|m|^3)/(1 + b1|m| + b2|m|^2 + b3|m|^3 + b4|m|^4)
	where |m| is the rotational running index, as defined previously above.

==============================================================================

How to use the broadening Python scripts?

The broadening Python files are labeled according to molecule type; the molecule the broadening file is labeled for should not be used on a different molecule.
For instance, the CO (Carbon Monoxide) broadening file should not be used to apply broadening to an SO2 (Sulfur Dioxide) line list.

To run these broadening Python files on your computer, make sure you have Python and the broadening Python files downloaded on your local machine.
Once the broadening Python files and the Input files are in the same folder on your machine, you can then run the Python scripts on the command line (example below).

Example of running the Python script of CO.py with Input file CO.par and Output file Output_CO.par with Anaconda's command prompt tool:
(First locate the directory where the downloaded Broadening Python files are located): cd filename (to go to the file named filename)
										       cd.. (to go back to the previous folder you were in)
										       Python CO.py (once you are in the folder with the Python script and input file)
										       CO.par (You will be asked to type the input file name)
										       Output_CO.par (You will be asked to type the output file name)
							(If the Python script works correctly then the last message is given below (in this case for CO.py))
							end for calculation: output "160.par + gamma_He + n_He + gamma_H2 + n_H2 + gamma_CO2 + n_CO2"
(This message means that the original HITRAN .par input file is given as part of the output, with extra columns (on the right-hand side) containing broadening of Helium, 
temperature dependence of Helium, broadening of Hydrogen, temperature dependence of Hydrogen, broadening of Carbon Dioxide and temperature dependence of Carbon Dioxide)

For questions related to using the Broadening Python files or about using HITRAN data please email info@hitran.org.

==============================================================================
