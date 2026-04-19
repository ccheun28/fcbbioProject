# prediction model of relapse vs. non-relapse based on clonal information and mutations in top 20 mutated genes

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt


# Load the combined clonal and mutation data
data = pd.read_csv('myclone-output_combined/myclone_results_combined_genes.tsv', sep='\t')
patient_gene_matrix = pd.read_csv('myclone-output_combined/myclone_patient-gene_matrix.tsv', sep='\t')
# Get top 20 mutated genes
patient_gene_matrix = patient_gene_matrix.drop(columns=['sample_id'])
gene_freq = patient_gene_matrix.sum(axis=0).sort_values(ascending=False)
top_genes = gene_freq.head(20).index

# Aggregate data per sample
aggregated = data.groupby('sample_id').agg({
    'cluster_id': 'nunique',  # clone count
    'shannon_index': 'first',  # shannon diversity (assuming same per sample)
    'MRD_EOI_Pct': 'first',  # MRD
    'relapse_status': 'first',  # relapse status
}).reset_index()
# include if any of the top 20 genes from top_genes are mutated in the sample separately
for gene in top_genes:
    aggregated[f'{gene}_mutated'] = patient_gene_matrix[gene].apply(lambda x: 1 if x > 0 else 0)

# Convert relapse_status to binary
aggregated['relapse'] = aggregated['relapse_status'].map({'Yes': 1, 'No': 0})

# Drop rows with missing values
aggregated = aggregated.dropna()

# Features and target
features = ['cluster_id', 'shannon_index', 'MRD_EOI_Pct'] + [f'{gene}_mutated' for gene in top_genes]
X = aggregated[features]
y = aggregated['relapse']

# Split data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

# Scale features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train logistic regression model
model = LogisticRegression(random_state=42)
model.fit(X_train_scaled, y_train)

# Predict
y_pred = model.predict(X_test_scaled)
y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]

# Evaluate
accuracy = accuracy_score(y_test, y_pred)
auc = roc_auc_score(y_test, y_pred_proba)
print(f'Accuracy: {accuracy:.3f}')
print(f'AUC: {auc:.3f}')
print('Classification Report:')
print(classification_report(y_test, y_pred))

# Feature importance (coefficients)
coefficients = model.coef_[0]
feature_importance = pd.DataFrame({
    'feature': features,
    'coefficient': coefficients
})
print('Feature Importance:')
print(feature_importance.sort_values('coefficient', ascending=False))

# Plot feature importance
plt.figure(figsize=(8, 6))
plt.barh(feature_importance['feature'], feature_importance['coefficient'])
plt.xlabel('Coefficient')
plt.title('Feature Importance in Logistic Regression')
plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/myclone_prediction_feature_importance.png')

# Plot ROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, label=f'ROC curve (area = {auc:.3f})')
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend(loc='lower right')
plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/myclone_prediction_roc_curve.png')