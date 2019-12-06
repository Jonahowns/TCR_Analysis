import numpy as np
import dcamethods as dca
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn import metrics
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn import preprocessing
from sklearn.metrics import confusion_matrix
import random


# DATA PREP
df = '/home/jonah/Dropbox (ASU)/Projects/DCA/ThrombinAptamers/v3/8all.txt'
affcutoff = 2

a, seqs = dca.Fasta_Read_Aff_wC(df, affcutoff)
X_pre, Y_pre = seqs, a

# Classify binders into two groups aff=2 and aff>2
def rewriteY_data(Y):
    newY = []
    for yid, y in enumerate(Y):
        if 1 < y < 3:
            newY.append(1.0)
        elif 3 <= y < 30:
            newY.append(2.0)
    return newY

# Count of Seqs in each group
Y_adj = rewriteY_data(Y_pre)
yaff = set(Y_adj)
for a in yaff:
    print('Aff:', a, 'Count:', Y_adj.count(a))

# One Hot Encoding Sequences
enc = preprocessing.OneHotEncoder()
X_all = []
for x in X_pre:
    X_all.append(list(x))
enc.fit(X_all)
onehotlabels = enc.transform(X_all).toarray()

# Separating our Training and Testing Data
X_train, X_test, Y_train, Y_test = train_test_split(onehotlabels, Y_adj, test_size=0.2, random_state=0)



# Random Forest Classifier
# gsc = GridSearchCV(estimator=RandomForestClassifier(),
#                    param_grid={'max_depth': range(8, 12), 'n_estimators': (1500, 2000, 2500, 3000),
#                                'max_features': (42, 44, 46)}, cv=5, scoring='neg_mean_squared_error', verbose=1, n_jobs=-1)
# grid_result = gsc.fit(X_train, Y_train)
# best_params = grid_result.best_params_
best_params = {'max_depth': 9, 'max_features': 46, 'n_estimators': 3000}
# print(best_params)

# rfr = RandomForestClassifier(max_depth=best_params["max_depth"], n_estimators=best_params["n_estimators"],
#                             max_features=best_params["max_features"], criterion="gini", random_state=False,
#                             verbose=True)

# scores = cross_val_score(rfr, X_train, Y_train, cv=10, scoring='neg_mean_absolute_error')
# rfr.fit(X_train, Y_train)
# ypred = rfr.predict(X_test)
# print(metrics.classification_report(ypred, Y_test))


# K Nearest Neighbors
kneb = KNeighborsClassifier(n_neighbors=5, weights='distance', algorithm='brute', p=2)
scores = cross_val_score(kneb, X_train, Y_train, cv=5, scoring='neg_mean_absolute_error')
print(scores)
kneb.fit(X_train, Y_train)
kpred = kneb.predict(X_test)
print(metrics.classification_report(kpred, Y_test))