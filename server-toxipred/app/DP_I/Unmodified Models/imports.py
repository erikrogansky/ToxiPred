# RDKit: https://www.rdkit.org/
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, MACCSkeys
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.ML.Descriptors import MoleculeDescriptors

# Pandas: https://pandas.pydata.org
import pandas as pd

# NumPy: https://numpy.org/doc/stable/release.html
import numpy as np

# Scikit-learn: https://scikit-learn.org
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold, RandomizedSearchCV
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import (accuracy_score, confusion_matrix, classification_report, 
                             roc_curve, roc_auc_score)
from sklearn.preprocessing import MinMaxScaler, FunctionTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.decomposition import KernelPCA
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import LabelEncoder
from sklearn.impute import SimpleImputer
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import precision_score, recall_score, f1_score

# XGBoost: https://xgboost.ai/
from sklearn.base import clone
import xgboost as xgb

# SHAP: https://github.com/slundberg/shap
import shap

# XSmiles: https://github.com/mahmoudnafifi/xsmiles
#import xsmiles

# JSON: https://docs.python.org/3/library/json.html
import json

# Matplotlib: https://matplotlib.org/
import matplotlib.pyplot as plt

# SciPy: https://www.scipy.org/
from scipy.stats import uniform

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import skew, kurtosis, shapiro
from sklearn.preprocessing import PowerTransformer, QuantileTransformer

import pandas as pd
from sklearn.decomposition import PCA

import optuna
import xgboost as xgb
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

from sklearn.linear_model import Lasso
import numpy as np

from sklearn.linear_model import LogisticRegression
import optuna
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, log_loss
from sklearn.metrics import accuracy_score, roc_curve, auc
from sklearn.model_selection import learning_curve


import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier


from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_curve, auc
from sklearn.base import clone
import matplotlib.pyplot as plt
import numpy as np

import joblib
import numpy as np
import pandas as pd
from collections import Counter
import torch
import torch.nn as nn

import xgboost as xgb
import optuna
from sklearn.model_selection import cross_val_score, StratifiedKFold, train_test_split
from sklearn.metrics import accuracy_scores