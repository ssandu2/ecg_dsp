# ecg_dsp
Collection of digital signal processing work on ECGs, ranging from segmentation to genetic fingerprinting. Work done at UICOM under Dr. Dawood Darbar and Dr. Anish S. Shah at UIC College of Medicine.
# Gene lists: 
1. AF_CM_genelist_R contains the script for classifying genes as sarcomeric, ion channel, or non-sarcomeric and non-ion channel (3 strata). 
2. Excel sheets contain 2 different gene lists, outputted from the script, that are labelled with classification and any notes, which are manually inputted.
3. NCBI Gene was used to download the list of genes found in AF_CM_Genelist1 and NCBI dbGaP was used to download the list of genes found in AF_CM_Genelist2 from Yoneda, Zachary T., et al. "Early-onset atrial fibrillation and the prevalence of rare variants in cardiomyopathy and arrhythmia genes." JAMA cardiology 6.12 (2021): 1371-1379.
4. The final combined gene list is AF_CM_Genelist_1+2_combined.
# F wave extraction:
1. Contains code for isolating ECGS in AF from LUDB dataset. The code to construct the ecgs, while modified, is largely from cardio_ml /LUDB_input_fft_fsst.
2. F wave extraction is done by two methods: FFT and PCA. In the frequency domains, f waves are extracted in the 4-10 hz band. PCA f wave extraction is in progress and aims to use PCA on raw ECG signal to output PCs corresponding to QRS complexes.
3. Goal is to wrap f wave extraction into a function with the argument of annotated or non-annotated ECG signal.
