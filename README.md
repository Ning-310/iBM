# Inference of Biomarker combinations with Minimal bias (iBM)

**Note: We developed a new computational pipeline named iBM, which consisted of three steps, including mutual DEPs or DEMs selection (MDS), candidate combination generation (CCG) to randomly select 10,000 biomarker combinations, and final combination prioritization (FCP) to get the protein or metabolite combination with a maximal accuracy and a minimal bias from the 5-fold cross-validation. The accuracy of a model was evaluated by calculating the total area under curve (AUC) value, and we also computed the total root mean squared error (RMSE) to measure the prediction bias. In the step of FCP, a widely used machine learning algorithm, penalized logistic regression (PLR), was used for model training and parameter optimization. The CC-specific biomarker combinations were separately determined for the proteomic and metabolic data.

## Requirements

The main requirements are listed below:

* Python 3.5
* Numpy
* Scikit-Learn
* Joblib
* Keras
* Pandas


## The description of iBM source codes

* DEPs or DEMs.py

    The code is used to select DEPs or DEMs.

* iBM.py

    The code is used to generate 10,000 groups of candidate biomarker combinations, and get the final biomarker combination.

* Train_biomarker.py

    The code is used to train the model for biomarker combinations.

* Test_biomarker.py

    The code is used to test the model for biomarker combinations.

* ROC.py

    The code is used to illustrate the receiver operating characteristic (ROC) curve based on sensitivity and 1-specificity scores, and compute the AUC value.


## The models in POC-19

* Protein_CC vs HC.model 

    The model is used for the classification of CC and HC by the protein biomarker combination.

* Protein_CC vs AC.model 

    The model is used for the classification of CC and AC by the protein biomarker combination.

* Metabolite_CC vs HC.model 

    The model is used for the classification of CC and HC by the metabolite biomarker combination.

* Metabolite_CC vs AC.model 

    The model is used for the classification of CC and AC by the metabolite biomarker combination.
