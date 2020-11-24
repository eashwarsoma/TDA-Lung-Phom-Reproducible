# TDA-Lung-Phom-Reproducible
Reproducible repository for "Persistent Homology of Tumor CT Scans Predicts Survival In Lung Cancer"

To reproduce results, first clone this repository to your computer. Next, download the Radiogenomics_DICOM and Radiomics_DICOM folders at these links (insert links). Move these folders to the same directory as the cloned repository. These DICOMs are directly from the TCIA and have been rearranged into a folder structure compatible with the code. 

To replicate the results, run the following 5 scripts in the stated order:
1. Segment_Script_1.R
2. Create_Numpy_Script_2.R
3. Homology_Make_Script_3.py
4. Format_Hom_Script_4.R
5. New_Analysis_Figures.R

For all of the scripts, you will likely have to install the packages and library before being able to run the code. All of the libraries and packages are named at the top of the script. 

If you have any questions regarding this repository, reproducing the results, or the paper, please don't hesitate to contact Eashwar Somasundaram at evs27@case.edu. 
