# Readme for planetary-broadeners

## License

See the LICENSE file available on this GitHub page.

These HITRAN Broadening Python files are distributed under the Apache 2.0 license

Copyright 2022 HITRAN Team. You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0


## Citation

Y. Tan, F.M. Skinner, S. Samuels, R.J. Hargreaves, R. Hashemi, I.E. Gordon (2022), Submitted to the Astrophysical Journal Supplement Series in March 2022
"H<sub>2</sub>, He, and CO<sub>2</sub> pressure-induced parameters for the HITRAN database. Part II: Line lists of CO<sub>2</sub>, N<sub>2</sub>O, CO, SO<sub>2</sub>, OH, OCS, H2CO, HCN, PH<sub>3</sub>, H<sub>2</sub>S and GeH<sub>4</sub>"


## Purpose

The broadening Python files are designed to apply broadening parameters to HITRAN line lists through dependency on rotational quantum numbers.


## Explanation of the Input & Output Files

Input files are given in the "Input-Broadening-Files" directory and are *sampled* portions of the HITRAN2020 line list data (not the complete line list).
In particular, the input files are HITRAN line lists in default format called .par (see [HITRAN Documentation](https://hitran.org/docs/definitions-and-units/) for an explanation of the .par format)

We highly recommend that users directly download the current line lists from the [HITRAN](https://hitran.org/) database prior to applying these broadening codes.

These broadening codes can be used on any line list, however the user must remember to alter the format of the read-in file if they are not using a HITRAN .par formatted line list.

	If the user is not using a HITRAN .par line list, then the primary columns to isolate when formatting 
	the read-in line list are the rotational quantum number(s) and the branch columns (for some molecules 
	the branches are not needed). This is due to the fact that the main data columns used to apply 
	broadening parameters are dependent on rotational quantum number (J lower, and for some molecules Ka 
    lower) and Branch letter (P, Q or R).

Output files are available in the "Output-Broadening-Files" directory. The output files are examples of what users should receive if they run the given input files correctly from the "Input Broadening Files" folder by using the broadening Python scripts, given in the "Broadening .py Files" folder.


## Downloading Broadening Parameters via HITRAN*online*

To access these additional foreign-broadening parameters via the HITRAN database directly, HITRAN users can proceed to the [HITRAN Database](https://hitran.org/) to create a customized output file when downloading line-by-line data.

Here is an example of adding H<sub>2</sub> broadening parameters to a line list output within the HITRAN database.
The steps are given as follows:

- Login (only requires your name and email address)
- Go to Line-by-Line under Data Access
- Choose a molecule
- Choose the isotopologues of the molecule
- Specify your spectral range
- (!STOP! This step is the important part where HITRAN is asking what the output format should be)
- Scroll down and click the green button (left-hand side) that says "Create New Output Format"
- (Now you should see a page that says "New Output Format")
- Under "Available Parameters" choose ".par line"
- Then under "Available Parameters" choose "&gamma;<sub>H2</sub>" then "n<sub>H2</sub>" and "&delta;<sub>H2</sub>"
- Give your output format a name within the "Output Format Name:" box
- Give your output format a description within the "Description:" box
- (You can also change the field separator and line endings, also you can apply a fixed width format and an output header line)
- Now click the green button "Save and Return to Data Search"
- (You might get a warning that reads: "You have not selected "Global isotopologue ID" or "Molecule ID" and "Isotopologue ID"  among of your parameters and so may have no way to know which transition corresponds to which species. Proceed?" Do not worry about this warning, we have the both the "Molecule ID" and "Isotopologue ID" within the HITRAN default .par format, which we already added to our output. Therefore, please proceed with confidence even if you receive this warning.)
- Press "Okay" to exit the warning
- Press "Start Data Search" to retrieve the final line list and its corresponding references with the broadening parameters included.


## Definitions

|m| values are the rotational running index. |m| is defined with the following relationships to the J rotational quanta:

- P-branch: m = -J"
- Q-branch: m = J"
- R-branch: m = J"+1

Broadening parameters refer to the following:

1. Lorentzian half-widths at half maximum (HWHM) denoted as &gamma;<sub>H2</sub>, &gamma;<sub>He</sub>, &gamma;<sub>CO2</sub> and &gamma;<sub>H2O</sub> for H<sub>2</sub>-, He-, CO<sub>2</sub>- and H<sub>2</sub>O-broadening respectively.
	
2. Temperature dependence (exponent) of these half widths, denoted as n<sub>H2</sub>, n<sub>He</sub>, n<sub>CO2</sub>, n<sub>H2O</sub> and defined through the power law (see below).
	
3. Collisional line shifts, denoted as &delta;<sub>H2</sub>, &delta;<sub>He</sub>, &delta;<sub>CO2</sub> and &delta;<sub>H2O</sub> which at the moment are available only for some HITRAN molecules and in some cases only for some specifically selected lines of these molecules. 

The power law equation for determining the HWHM at T (Temperature) is given as: &gamma;(T) = &gamma;(T<sub>0</sub>)([T<sub>0</sub>/T])<sup>n</sup> where T<sub>0</sub> is the reference temperature (296K in HITRAN) and &gamma;(T<sub>0</sub>) is the HWHM at the reference temperature.

The broadening Python files utilize 3<sup>rd</sup>-to-4<sup>th</sup> order Pad&eacute; approximants (equation given below) to populate the broadening data throughout the line list.

The Padé approximants were fit to available laboratory data sets in order to extrapolate to other transitions where data was not available.

The 3<sup>rd</sup>-to-4<sup>th</sup> order Pad&eacute; approximant: 
- &gamma;<sub>x</sub>(|m|) = (a<sub>0</sub>+a<sub>1</sub>|m|+a<sub>2</sub>|m|<sup>2</sup>+a<sub>3</sub>|m|<sup>3</sup>)/(1+b<sub>1</sub>|m|+b<sub>2</sub>|m|<sup>2</sup>+b<sub>3</sub>|m|<sup>3</sup>+b<sub>4</sub>|m|<sup>4</sup>) 
where |m| is the rotational running index, as defined previously above.


## How to use the broadening Python scripts?

The broadening Python files are labeled according to molecule type; the molecule the broadening file is labeled for should not be used on a different molecule.
For instance, the CO (Carbon Monoxide) broadening file should not be used to apply broadening to an SO$_2$ (Sulfur Dioxide) line list.

To run these broadening Python files on your computer, make sure you have Python installed and the broadening Python files downloaded on your local machine.

- Download the complete "planetary-broadeners" repository directory to local storage from GitHub
- Alternatively, clone the "planetary-broadeners" git repository
        
Once the broadening Python files and the Input files are downloaded, you can then run the Python scripts on the command line (example below).
```
# Example of running the Python script CO.py with input file sample_CO.par to create the output  
# file sample_CO_out.par from the command line:
	
cd /full-path/Broadening-Files              # Move to the Broadening-Files directory, in the
                                            # "full-path" on your local machine

python CO.py                                # Run the Python script for carbon monoxide

Input-Broadening-Files/sample_CO.par        # You will be asked to enter the input filename. 
                                            # This example uses the sample_CO.par file in the  
                                            # Input-Broadening-Files directory

sample_CO_out.par                           # You will be asked to enter the output filename. 
                                            # This example uses recreates the sample_CO_out.par   
                                            # file from the Output-Broadening-Files directory
							       
# If the Python script is successful, then a final output message will be given (e.g., the 
# message below is from CO.py)

> end of calculation: output "160.par + gamma_He + n_He + gamma_H2 + n_H2 + gamma_CO2 + n_CO2"

# This message means that the original HITRAN .par input file is given as part of the output 
# (160.par), with additional columns containing the pressure broadening due to helium (gamma_He), 
# temperature dependence of helium broadening (n_He), broadening due to hydrogen (gamma_H2), 
# temperature dependence of hydrogen broadening (n_H2), broadening due to carbon dioxide 
# (gamma_CO2), and temperature dependence of carbon dioxide broadening (n_CO2).
```
For questions related to using the Broadening Python files or about using HITRAN data please email info@hitran.org.



