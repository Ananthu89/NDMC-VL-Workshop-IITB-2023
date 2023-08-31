# General instructions 

Here you will find a package in `R` for modelling the visceral leishmaniasis transmission cycle in India. 

1. **Read this first and save your time!.** <br/> 
You can see five Sections (I-V) below. Section I is for everyone. Section II is not relevant for those familiar with `RStudio` and `R` and have already installed the `pomp` package (as recommended in advance). Section III is for those unfamiliar with `R` and `C`. Sections IV and V are meant for those who want to adapt the scripts for their specific needs.  

2. **Meant for whom?** <br/> 
Even for those without any understanding of programming languages, if they have a quantitative background (a masters or 4-year bachelors) and/or basic research experience in mathematical modelling, it should be possible to follow the instructions and use the scripts in the existing forms or adapt for their needs. 

3. **Some (familiar) advice.** <br/> 
Before making changes in any of the scripts, we strongly advise you to save each file as a separate copy. In that case, if at all you encounter an error due to the modification you may make, you can compare with the original versions you received from us.  


-------------------------------------------------------------------

## I. Files & sub-folders and their arrangements 


NOTE: Run the scripts in the order from 01 to (scripts in folder) 04. There are two reasons for that. 

1. This (01->04) is the order of increasing complexity, as you can also see from the Table below, which provides details on the purpose of each script and how they are different from each other.   

2. Script 02 needs to be run to generate the two folders where the output figures will be stored. Script 03 and the scripts in folder 04 will run only if those folders are already present. (For these scripts, there will be two curves per the output plot, red line representing deterministic output and blue curve showing stochastic output.) 


| Sub-folder    | File    | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Contents&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Additional&nbsp;comments&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     |
|:-------------:|:-------:|:---------:|:-----------------------:|
| 00_Source     | param_values.R | Script containing parameter values | Directly called in the scripts in the main folder or other sub-folders. Hence, no need to run this separately. |
| 00_Source     | model.R        | Script containing deterministic and stochastic models | Directly called in the scripts in the main folder or other sub-folders. Hence, no need to run this separately. |
| -      | 01_parameter_checks_example.R | Script to perform some basic tests before running the full models | Expected values for most of the tests are given there. |
| -     | 02_fixed_foi.R        | Script to run **_model.R_**, using the values in **_param_values.R_**, and process and plot the output, assuming the infection is sustained only through incoming sandflies from neighbourhood (known as the external force of infection) | Transmission-dynamics of sandflies is absent here. The incoming flies can transmit infection to humans, but humans cannot infect them back. (This script has two versions: **population-model** (total population size being constant) and **cohort-model** (total population size decreasing over time).) | 
| -      | 03_sandfly_tran_dyn.R | Script to run **_model.R_**, using the values in **_param_values.R_**, and process and plot the output, assuming transmission between humans and sandflies (sandfly-transmission dynamics) is present | The sandfly-abundance follows a seasonal pattern. External force of infection is not needed to sustain transmission, and it is assumed to be absent. (When sandfly-transmission is present, the model used is always the **population-model**.) |
| 04_Interventions      | policy_highdetection.R | **_03_sandfly_tran_dyn.R_**, with the parameter pertaining to improved detection of symptomatic VL cases changed at a certain time-point, thus representing a policy/intervention to control the infection | Use of *covar* argument in `pomp` |
| 04_Interventions      | policy_fly_abundance_reduction.R | **_03_sandfly_tran_dyn.R_**, with the parameter pertaining to reduction in sandfly abundance (through efforts such as IRS (indoor residual spraying)) changed at a certain time-point, thus representing a policy/intervention to control the infection | Use of *covar* argument in `pomp` |
| -      | Model description.pdf | Document containing schematic diagram, equations, and details regarding the parameter values including references | This file is NOT required to run any of the scripts. It's for the readers' reference. |
| -      | VL modelling workshop IIT - Aug 2023_Day1.pdf | The slides that have been presented on day 1 | This file is NOT required to run any of the scripts. It's for the readers' reference. |
| -      | VL modelling workshop IIT - Aug 2023_Day2.pdf | The slides that have been presented on day 2 | This file is NOT required to run any of the scripts. It's for the readers' reference. |

-------------------------------------------------------------------

## II. Requirements to save and run the scripts 


This section (except step 2) is particularly for those who are unfamiliar with RStudio or R. Others may skip this. 


**Saving the files** 


Save them in the same way they appear in the folders here. Also, to read the instructions presented here, open this **README.md** file on Github and not on your system. 


**Running the scripts**


1. `RStudio` and `R`  

2. `pomp` package (see the instructions for the installation: <https://kingaa.github.io/sbied/prep/>) 

3. R-libraries `ggplot2`, `rstudioapi`, `gridExtra`, and `data.table`  

The below commands will be useful in case you need to install the required R-libraries to run the scripts:

_install.packages("ggplot2")_
_install.packages("rstudioapi")_   
 _install.packages("data.table")_
_install.packages("gridExtra")_

4. Finally, to run a script, after opening it using `RStudio`, select the portion of the script you want to run and click on the **`Run`** button (at the right top) of the screen 

------------------------------------------------------------------- 

## III. Languages used and 'comments' in the scripts   

1. All the scripts except **_model.R_** is written exclusively using `R`. In a line, any text that is to the right of the symbol **"#"** is a comment in R, meant only for the reader's understanding. Such portions are always excluded when we run a script in R. 

2. The script **_model.R_** is written in `R`, but four of the R-objects (*vl_code_det*, *vl_det_init*, *vl_code_sto*, *vl_sto_init*) have portions written in C (known as C-Snippets) inside inverted commas (" "). In a line, any text that is to the right of **"//"** is a comment, meant only for the reader's understanding. These are also excluded when the script is run. 

-------------------------------------------------------------------

## IV. When a transition is added (removed) between two compartments 


A transition can be added (removed) with or without adding (removing) a parameter.  


### Addition (removal) of a transition without changing parameters  

This is a relatively easy process. Only the script **_model.R_** needs to be adapted in this case. The specific details of the required changes are listed below. The listed steps are also applicable when a parameter change is involved. 

Object&nbsp;in&nbsp;R&nbsp;to&nbsp;be&nbsp;changed   | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Relevant&nbsp;details&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; |  
|:-----------:|:-----------------------:| 
| vl_code_det | The differential equations involving the added (removed) transition. (The differential equation for S might also have to be changed.) | 
| vl_code_sto | 1. For the compartment 'Z' to (from) which the transition needs to be added (removed), the dimensionality 'n' of both the concerned transition-number (representing the number of individuals exiting a compartment), denoted as *e_Z[n]*, and transition-rate, denoted as *rate_Z[n]*, increases (decreases) by *1*. <br/> 2. Definition of the added (removed) transition-rate. A dimensionality-change needs to be made in <br/> *reulermultinom(n, Z, &rate_Z[0], dt, &e_Z[0]);*. <br/> 3. The system of ODEs involving the added (removed) transition-numbers. (This may also involve the accumulators for the VL incidences.) | 


### Addition (removal) of a transition involving a parameter addition (removal)    


1. In this case, apart from the above listed aspects, the change needs to be implemented in the R-object "param_values" of **_param_values.R_**. 

2. Other changes required in the **main scripts** and **_01_parameter_checks_example.R_** are dependent on whether this parameter is needed (used) in the these scripts. This list is too long to be specifically pointed out. Some useful suggestions are given below. <br/> (1) If the added (removed) parameter *w* needs to be (is) changed in the **main scripts** using the command *param_values["w"]*, a change is required here. <br/> (2) In the case of removal of a parameter, search for this parameter in the scripts and modify this part accordingly. <br/> (3) For addition of a parameter, changes are applicable if those parameters affect the VL-incidence calculation and if those need to be changed as part of an intervention.  


NOTE: The addition (removal) of a transition is not necessarily accompanied by the addition (removal) of a compartment. On the other hand, the addition (removal) of a compartment is always accompanied by that of one or more transitions. 


-------------------------------------------------------------------

## V. Checklist for the addition or removal of a compartment 


The addition (removal) of a compartment affects all the scripts. The below Table presents a documentation of all the possible changes you need to make when you change a compartment. All these changes are not always required, as mentioned in the third column. The last column provides useful info such as the kind of change required and the reason, with repetitions of this change being highlighted. (A **main script** refers to **_02_fixed_foi.R_**, **_03_sandfly_tran_dyn.R_**, or the scripts in the folder *`04_Interventions`*.) 


| Script    |  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Object&nbsp;or&nbsp;function&nbsp;in&nbsp;R&nbsp;to&nbsp;be&nbsp;changed&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  <br/> (exact argument/command given as it is) | Always needed? | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Relevant&nbsp;details&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | 
|:-------------:|:----------------:|:----:|:-----------------------:| 
| model.R    |  *statenames_det*   | Yes | Addition (removal) of the compartment name | 
| model.R    |  *vl_code_det*   | Yes | The expression for *N* and the differential equations involving the added (removed) compartment. (The expression for *lambdaHF* and the differential equation for *S* might also have to be changed.) | 
| model.R    |  *vl_det_init*   | Yes | The expressions for *N_0* and the added (removed) compartment | 
| model.R    |  *vl_code_sto*   | Yes | 1. Declarations of both the transition-number (representing the number of individuals exiting a compartment), denoted as *e_Z[n]*, and transition-rate, denoted as *rate_Z[n]*, involving the added (removed) compartment 'Z', with 'n' representing the dimensionality. <br/> 2. Definitions of *N* and maybe *lambdaHF* as well. <br/> 3. Definitions of the transition-terms for the added (removed) compartment. <br/> 4. The system of ODEs involving the transition-numbers for the added (removed) compartment. (This may also involve the accumulators for the VL incidences.) | 
| model.R    |  *vl_sto_init*   | Yes | The expression for the added (removed) compartment | 
| param_values.R    |  *init_values*   | Yes | Initial value for the added (removed) compartment |
| Main script |  *names(output_det)[1:13]*   | Yes | Reason: total number of compartments would have increased (decreased) from 13. (Note that this may appear **twice** if VL-incidence needs to be plotted as a function of sandfly abundance.) |
| Main script |  *output_det[, I := I11 + I12 + I13 + I21 + I22 + I23]*   | No | Only needed if the added (removed) compartment is part of the VL-symptomatic ones. |
| Main script |  *output_det[, Ixy_detected := rhoIx * Ixy]*, <br/> where *x=1,2* and *y=1,2,3*  | No | Only needed if the added (removed) compartment is part of the VL-symptomatic ones (*Ixy*). (Currently, there are six such lines for different combinations of *x, y*. One line might have to be added or removed. Note also that this may appear **twice** if VL-incidence needs to be plotted as a function of sandfly abundance.) |
| Main script |  *output_det[, VL_symp_inc := (fs * rhoH * H) / N]*  | No | Needs to be removed if the removed compartment is H. (Warning: this is an unlikely scenario, as the removal of this compartment changes the system completely.)  |
| Main script |  *output_det[, VL_death_1 := 3 * muVL * I13]*   | No | A change is needed here only if the addition (removal) of compartment affects this part directly. |
| Main script |  *output_det[, VL_death_2 := 3 * muVL * I23]*   | No | A change is needed here only if the addition (removal) of compartment affects this part directly. |
| Main script |  *names(output_sto)[3:15]*   | Yes | Reason: total number of compartments would have increased (decreased) from 13. |
| Main script |  *output_sto[, I := I11 + I12 + I13 + I21 + I22 + I23]*   | No | Only needed if the added (removed) compartment is part of the VL-symptomatic ones. |
| Main script |  *plot_Z*   | No | Only needed if the added (removed) compartment *Z* is being plotted. |
| Main script |  *ggsave(...)* | No | Only needed if the plot for the added (removed) compartment is being saved. |

  