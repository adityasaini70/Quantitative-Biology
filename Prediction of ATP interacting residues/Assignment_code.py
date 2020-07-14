#!/usr/bin/env python
# coding: utf-8

# # IQB Assignment - 3 : Prediction of ATP interacting residues
# ### Aditya Saini, 2018125, ECE
# ---
# 
# ### Objective
# The objective of this assignment was to detect the presence of ATP interacting residues in a given protein sequence. **The main references of this assignment are the videos uploaded by Sir on Google classroom and** [**this research paper**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-434).
# 
# 
# ### Directory structure
# The structure of the project directory is as follows:
# 
# * **Assignment_notebook.ipynb :** Jupyter notebook implementation of the Assignment
# * **Assignment_code.py :** Python code for the assignment
# * **Report.pdf :** A pdf version of _Assignment_notebook.ipynb_. Contains information regarding usage of code, implementation details, etc.
# * **train.data :** Data used for creation of dataset
# * **test1.txt :** Test data
# * **out_new.csv :** Output data used for submission on Kaggle
# 
# ### Usage
# Any Python 3.x version will work for the implementation. 
# > For running the script:
# `
#     python Assignment_code.py
# `
# 
# ### Workflow
# The workflow for this assignment is as follows:
# 1. [**Sec-1 : Extracting data out of 'train.data' to create a viable dataset**](#Sec-1-:-Extracting-data-out-of-__train.data__-to-create-a-viable-dataset)
# 1. [**Sec-2 : Creating a dataset from extracted data**](#Sec-2-:-Creating-a-dataset-from-extracted-data)
#     1. [**Sec-2.1 : Generating patterns of given sequences**](#Sec-2.1-:-Generating-patterns-of-given-sequences) 
#     1. [**Sec-2.2 : Creating labels for our data**](#Sec-2.2-:-Creating-labels-for-our-data) 
#     1. [**Sec-2.3 : Creating a binary profile for our generated patterns**](#Sec-2.3-:-Creating-a-binary-profile-for-our-generated-patterns)
# 1. [**Sec-3 : Applying Machine Learning techniques to dataset**](#Sec-3-:-Applying-Machine-Learning-techniques-to-dataset)
#     1. [**Sec-3.1 : Splitting into training and testing sets**](#Sec-3.1-:-Splitting-into-training-and-testing-sets) 
#     1. [**Sec-3.2 : Fitting our data on the training sets**](#Sec-3.2-:-Fitting-our-data-on-the-training-sets) 
#     1. [**Sec-3.3 : Predicting the test data**](#Sec-3.3-:-Predicting-the-test-data)
#     1. [**Sec-3.4 : Evaluating our model by calculating the area under ROC curve**](#Sec-3.4-:-Evaluating-our-model-by-calculating-the-area-under-ROC-curve)
# 1. [**Sec-4 : Generating output for submission**](#Sec-4-:-Generating-output-for-submission)
# 

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm.auto import tqdm


# # Sec-1 : Extracting data out of __train.data__ to create a viable dataset

# In[ ]:


data = pd.read_csv('train.data')


# In[3]:


data


# # Sec-2 : Creating a dataset from extracted data 

# In[ ]:


sequence_data = data['Amino Acid Sequence']


# In[5]:


sequence_data


# # Sec-2.1 : Generating patterns of given sequences 

# In[ ]:


def generate_pattern(seq_data, window_size):
    dummy_variable_length = int((window_size-1)/2)
    ans = [] #Will hold the different patterns for a given amino acid sequence
    
    for sequence in seq_data:
        #Adding dummy variables to the extreme ends of the string
        string = "X" * dummy_variable_length + sequence + "X" * dummy_variable_length

        #Generating the patterns of size = "window_size"

        for idx in range(0, len(string) + 1 - window_size):
            ans.append(string[idx : idx + window_size])
            
    return pd.Series(ans)


# #### Referring the research paper and from cross-validation, I set the window-size to 17, for achieving the best results.
# 

# In[ ]:


size = 17
pattern_data = generate_pattern(sequence_data, size)


# In[10]:


pattern_data


# # Sec-2.2 : Creating labels for our data
# 
# * If the middle element was **lower case**, a **+1** was assigned to the sequence, denoting that the given **residue is an ATP interacting residue**
# * If the middle element was **upper case**, a **-1** was assigned to the sequence, denoting that the given **residue is NOT an ATP interacting residue**

# In[ ]:


def label_data(pat_data):
    ans = []
    for protein_seq in pat_data:
        target = protein_seq[int(len(protein_seq)/2)]
        if target.islower():
            ans.append(1)
        elif target.isupper():
            ans.append(-1)
    return pd.Series(ans)


# In[ ]:


Y_data = label_data(pattern_data)


# In[13]:


Y_data


# # Sec-2.3 : Creating a binary profile for our generated patterns
# 
# Each pattern sequence will be matched to a vector sequence, consisting of amino acids and, a 17\*21 length vector will be generated. This sequence is a _**reshaped (or flattened)**_ binary matrix which will help us in representing our patterns quantitatively.
# 

# In[ ]:


def generate_binaryProfile(pat_data):
    amino_acid = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X' ]    
    ans = []
    for series in tqdm(pat_data):
        ans.append([])
        for acid_1 in series:
            for acid_2 in amino_acid:
                if(acid_1.upper() == acid_2):
                    ans[-1].append(1)
                else:
                    ans[-1].append(0)
    
    return pd.Series(ans)


# In[15]:


X_data = generate_binaryProfile(pattern_data)


# In[16]:


X_data


# # Sec-3 : Applying Machine Learning techniques to dataset
# 
# We chose SVC (Support Vector Classifier) as our training algorithm as it is best suited for binary class classification problem on huge datasets.

# In[ ]:


from sklearn.svm import SVC
model = SVC(kernel = 'rbf', C = 2, gamma = 0.1)


# # Sec-3.1 : Splitting into training and testing sets
# 

# In[ ]:


from sklearn.model_selection import train_test_split
X_train, X_test, Y_train, Y_test = train_test_split(X_data.tolist(), Y_data.tolist(), test_size=0.2, random_state=0)


# # Sec-3.2 : Fitting our data on the training sets

# In[19]:


model.fit(X_train, Y_train)


# # Sec-3.3 : Predicting the test data

# In[ ]:


Y_test_predict = model.predict(X_test) 


# # Sec-3.4 : Evaluating our model by calculating the area under ROC curve

# In[21]:


from sklearn.metrics import roc_auc_score
roc_auc_score(Y_test, Y_test_predict.tolist())


# In[23]:


Y_test_predict


# # Sec-4 : Generating output for submission
# 

# In[ ]:


pd.read_csv('test1.txt')
X_predict_data = pd.read_csv('test1.txt')['Lable'].tolist()
X_predict_string = ''.join(map(str, X_predict_data))


# In[ ]:


X_predict_pattern_data = generate_pattern([X_predict_string], size)


# In[29]:


X_predict_pattern_data


# In[30]:


X_predict = generate_binaryProfile(X_predict_pattern_data)


# In[ ]:


Y_predict = model.predict(X_predict.tolist())


# In[32]:


pd.Series(Y_predict)


# In[33]:


output = { 'ID': pd.read_csv('test1.txt')['ID'], 'Lable': pd.Series(Y_predict) } 
output = pd.DataFrame(output)
output


# In[ ]:


output.to_csv('out_new.csv', index=False)  

