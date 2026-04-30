# prediction model of relapse vs. non-relapse based on clonal information and mutations in top 20 mutated genes

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

from sklearn.metrics import accuracy_score, classification_report, roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

myclone = 0 # set to 1 to run with myclone data, 0 for pyclone data
if myclone:
    data = pd.read_csv('myclone-output_combined/myclone_results_combined_genes.tsv', sep='\t')
    patient_gene_matrix = pd.read_csv('myclone-output_combined/myclone_patient-gene_matrix.tsv', sep='\t')
    fig_features = 'myclone_prediction_feature_importance_rf.png'
    fig_roc = 'myclone_prediction_roc_curve_new.png'
else:
    data = pd.read_csv('pyclone_output_combined/pyclone_output_combined_genes.tsv', sep='\t')
    patient_gene_matrix = pd.read_csv('pyclone_output_combined/pyclone_patient-gene_matrix.tsv', sep='\t')
    fig_features = 'pyclone_prediction_feature_importance_rf.png'
    fig_roc = 'pyclone_prediction_roc_curve_new.png'

# Load the combined clonal and mutation data
# Get top 20 mutated genes
patient_gene_matrix = patient_gene_matrix.drop(columns=['sample_id'])
# gene_freq = patient_gene_matrix.sum(axis=0).sort_values(ascending=False)
# top_genes = gene_freq.head(20).index
treatment_data = pd.read_csv('treatment_data.csv')

# gene sets
top_genes = ['KRAS'
    ,'NRAS'
    ,'CREBBP'
    ,'PAX5'
    ,'ETV6'
    ,'IKZF1'
    ,'TP53'
    ,'PTPN11'
    ,'FLT3'
    ,'JAK2'
    ,'NSD2'
    ,'KMT2D'
    ,'SETD2'
]
# Ras pathway genes
ras_genes = ['KRAS', 'NRAS', 'PTPN11', 'FLT3']
# driver genes for specific leukemia subtypes
driver_genes = ['ETV6', 'IKZF1', 'PAX5', 'JAK2', 'NSD2']


# top_genes minus ras_genes and driver genes
final_genes = [gene for gene in top_genes if gene not in ras_genes and gene not in driver_genes]

# Aggregate data per sample
aggregated = data.groupby('sample_id').agg({
    'cluster_id': 'nunique',  # clone count
    'shannon_index': 'first',  # shannon diversity (assuming same per sample)
    'MRD_EOI_Pct': 'first',  # MRD
    'relapse_status': 'first',  # relapse status
}).reset_index()
# add treatment info
aggregated = aggregated.merge(treatment_data, left_on='sample_id', right_on='sample_id', how='left')

# include if any of the top 20 genes from top_genes are mutated in the sample separately
for gene in final_genes:
    aggregated[f'{gene}_mutated'] = patient_gene_matrix[gene].apply(lambda x: 1 if x > 0 else 0)

# if any of the ras pathway genes are mutated, create a binary feature for that
aggregated['ras_pathway_mutated'] = patient_gene_matrix[ras_genes].apply(lambda x: 1 if x.sum() > 0 else 0, axis=1)
# if any of the driver pathway genes are mutated, create a binary feature for that
aggregated['driver_pathway_mutated'] = patient_gene_matrix[driver_genes].apply(lambda x: 1 if x.sum() > 0 else 0, axis=1)

# Convert relapse_status to binary
aggregated['relapse'] = aggregated['relapse_status'].map({'Yes': 1, 'No': 0})

# Drop rows with missing values
aggregated = aggregated.dropna()


# PREDICTION MODEL
print(aggregated.head())
# Features and target
features = ['cluster_id', 'shannon_index', 'MRD_EOI_Pct', 'ras_pathway_mutated', 'driver_pathway_mutated', 'risk_group', 'intensity', 'down_syndrome'] + [f'{gene}_mutated' for gene in final_genes]
X = aggregated[features]
y = aggregated['relapse']

# Split data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=26, stratify=y)

# Scale features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train random forest classifier
model = RandomForestClassifier(n_estimators=100, random_state=42)
# model = LogisticRegression(random_state=42)
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

# Feature importance (importance scores)
importance_scores = model.feature_importances_
feature_importance = pd.DataFrame({
    'feature': features,
    'importance': importance_scores
})
print('Feature Importance:')
print(feature_importance.sort_values('importance', ascending=False))

# Plot feature importance
plt.figure(figsize=(16, 10))
plt.barh(feature_importance['feature'], feature_importance['importance'])
plt.xlabel('Importance')
plt.title('Feature Importance in Random Forest Model')
plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/'+fig_features)

# Plot ROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
plt.figure(figsize=(12, 10))
plt.plot(fpr, tpr, label=f'ROC curve (area = {auc:.3f})')
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend(loc='lower right')
plt.savefig('/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/fcbbioProject/figures/'+fig_roc)