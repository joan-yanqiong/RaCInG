# RaCInG
The random graph model to infer cell-cell communication networks based on bulk RNA-seq data. The ropo contains the code used for the paper "Mathematically mapping the network of cells in the tumor microenvironment" by van Santvoort et al.

## General description
In this study we used patient specific bulk RNA-seq input together with non-patient specific prior knowledge on possible ligand-receptor interactions to reconstruct cell-cell interaction networks in an indivudal patient's tumour. The model operates in four main steps:
1. It transforms mRNA-seq input data into four matrices used by the graph generation procedure. These four input matrices contain information about cellular deconvolution (C-matrix), ligand-receptor interaction pair quantification (LR-matrix), ligand compatibility with cell-types (L-matrix) and receptor compatibility with cell types (R-matrix).
2. It uses the matrices to construct an ensemble of possible cell-cell interaction networks for individual patients. It does this by generating a fixed number of individual cells from the C-matrix, and a fixed number of ligand-receptor pair from the LR-matrix. Then, it binds the LR-pairs to indivudual cells unfiromly at random, respecting the compatibility of the ligand/receptor with the cell-type until everything is paired.
3. From the ensemble of networks certain fingerprints are extracted (e.g. the count of triangles). By averaging over the the counts for different graphs in the ensemble a feature value for the individual patient is created. Features are only extracted if they are expected to remain statistically consistent over the networks ensemble.
4. The features are used as biomarkers for the individual patient's tumor micro-environment, and through statistical testing meta-features of a patient (like resonse to immunotherapy) can be analysed. 

To validate the model, we have applied it to extract 444 network features related to the tumor microenvironment in 3213 solid cancer patients, and unveiled associations with immune response and subtypes, and identified cancer-specific differences in inter-cellular signaling. Additionally, we have used RaCInG in 118 patients with known response to immunotherapy to explain how immune phenotypes regulated by context-specific intercellular communication affect their response. RaCInG is a modular pipeline, and we envision its application for cell-cell interaction reconstruction in different contexts.

## Uses of the model
The code can be used to transform RNA-seq data into input matrices for RaCInG, and execute it to construct cell-cell interaction networks for individual patients. Moreover, it can be used to extract some predefined features for these networks, and do some statistical tests on these features. In order to run this model, the folder *R_code* contains all functionalities to turn RNA-seq data into input matrices for the model. The folder *Python_Code* contains the functions to execute the additional functionalities:
- The function *execute_sim.py* allows the user to construct networks based on the input matrices and extract certain features from them.
- The function *txt_to_csv.py* combines the output of a uniform and non-uniform run of the random graph model to create normalized feature values for each patient.
- The function *statistical_analysis.py* allows the user to do some statistical analysis on the output of the random graph model.
- The function *retrieveLigRecInfo.py* allows the user to compute the ligand-receptor interaction probabilities in given cell-types.
- The function *Circos.py* allows the user to compute direct communication values using the theoretical kernel of the random graph model to be plotted in a Circos plot.

In execution of the Python code make sure that all Python functions are located together. They each have their own functionalities, but work together to execute the model. Moreover, make sure that all input matrices for the model are located in the same folder (in our case the folder *Example input*) and that the metadata is located in the same folder as the *statistical_analysis.py* function.

<mark>A demo of the Python functionalities is available as Jupyter notebook in the Python_Code folder (called *Demo.ipynb*). A demo of R functionalities is available in the R_Code folder as Rmarkdown files.</mark>
