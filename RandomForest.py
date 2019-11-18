import numpy as np
import dcamethods as dca
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn import metrics
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn import preprocessing
from sklearn.metrics import confusion_matrix
import random

df = '/home/jonah/Dropbox (ASU)/Projects/DCA/ThrombinAptamers/v3/8all.txt'
affcutoff = 2

a, seqs = dca.Fasta_Read_Aff_wC(df, affcutoff)
X_pre, Y_pre = seqs, a

def rewrite_data(X, Y):
    for yid, y in enumerate(Y):


enc = preprocessing.OneHotEncoder()
X_all = []
for x in X_pre:
    X_all.append(list(x))
enc.fit(X_all)
onehotlabels = enc.transform(X_all).toarray()
print(onehotlabels[0])
print(X_all[0])


print(len(onehotlabels))
print(len(Y_pre))
X_train, X_test, Y_train, Y_test = train_test_split(onehotlabels, Y_pre, test_size=0.2, random_state=0)

print(Y_test)
# gsc = GridSearchCV(estimator=RandomForestClassifier(),
#                    param_grid={'max_depth': range(11, 18), 'n_estimators': (100, 250, 500, 750, 1000, 1500),
#                                'max_features': (8, 12, 18, 24)}, cv=5, scoring= 'neg_mean_squared_error', verbose=1, n_jobs=-1)
# grid_result = gsc.fit(X_train, Y_train)
# best_params = grid_result.best_params_
best_params = {'max_depth': 16, 'max_features': 18, 'n_estimators': 250}
# print(best_params)

rfr = RandomForestClassifier(max_depth=best_params["max_depth"], n_estimators=best_params["n_estimators"],
                            max_features=best_params["max_features"], criterion="gini", random_state=False,
                            verbose=True)

scores = cross_val_score(rfr, X_train, Y_train, cv=10, scoring='neg_mean_absolute_error')
print(scores)

rfr.fit(X_train, Y_train)
ypred = rfr.predict(X_test)
print(ypred)
print(Y_test)
print(metrics.classification_report(ypred, Y_test))
# mat = confusion_matrix(Y_test, ypred)
# sns.heatmap(mat.T, square=True, annot=True, fmt='d', cbar=False)
# plt.xlabel("true label")
# plt.ylabel("predicted label")
# plt.savefig('/home/jonah/Desktop/rftrial.png', dpi=600)








