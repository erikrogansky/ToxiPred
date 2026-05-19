export interface ModelMetric {
  label: string;
  train?: string;
  test?: string;
  stacking?: string;
  voting?: string;
}

export interface ModelDocumentation {
  modelName: string;
  slug: string;
  name: string;
  badge: string;
  endpoint: string;
  assayType: string;
  algorithm: string;
  dataset: string;
  observations: number;
  classBalance: string[];
  descriptors: string;
  fingerprint: string;
  featureSelection: string;
  tuning: string;
  threshold: string;
  homeDetails: string[];
  metrics: ModelMetric[];
  hyperparameters: string[];
  note?: string;
  documentationHref: string;
}

const docsBase = '/documentation#';

export const modelDocumentationByName: Record<string, ModelDocumentation> = {
  'XGB Corrosion': {
    modelName: 'XGB Corrosion',
    slug: 'model-xgb-corrosion',
    name: 'Corrosion 3D - In Vitro',
    badge: 'In Vitro',
    endpoint: 'Corrosion',
    assayType: 'In vitro 3D',
    algorithm: 'XGBoost classifier',
    dataset: 'In vitro 3D corrosion dataset',
    observations: 151,
    classBalance: ['Positive (1): 60 / 0.40', 'Negative (0): 91 / 0.60'],
    descriptors: '23 molecular descriptors',
    fingerprint: 'None',
    featureSelection: 'XGB Importance',
    tuning: 'Randomized Search + Manual Fine-Tuning',
    threshold: '0.5',
    homeDetails: [
      'XGBoost classifier',
      '23 molecular descriptors',
      '151-compound in vitro 3D corrosion dataset',
      'Test AUC 0.83, accuracy 0.81',
    ],
    metrics: [
      { label: 'Accuracy', train: '0.82', test: '0.81' },
      { label: 'AUC', train: '0.85', test: '0.83' },
      { label: 'F1', train: '-', test: '0.75' },
      { label: 'Recall', train: '-', test: '0.75' },
      { label: 'Precision', train: '-', test: '0.75' },
    ],
    hyperparameters: [
      'n_estimators: 30',
      'max_depth: 4',
      'learning_rate: 0.1',
      'subsample: 0.9',
      'colsample_bytree: 0.6',
      'min_child_weight: 1',
      'gamma: 0',
      'reg_lambda: 1',
      'reg_alpha: 0',
    ],
    documentationHref: `${docsBase}model-xgb-corrosion`,
  },
  'GB Corrosion (in vivo)': {
    modelName: 'GB Corrosion (in vivo)',
    slug: 'model-gb-corrosion-in-vivo',
    name: 'Corrosion - In Vivo',
    badge: 'In Vivo',
    endpoint: 'Corrosion',
    assayType: 'In vivo',
    algorithm: 'Gradient Boosting classifier',
    dataset: 'In vivo corrosion dataset',
    observations: 189,
    classBalance: ['Positive (1): 57 / 0.30', 'Negative (0): 132 / 0.70'],
    descriptors: '6 molecular descriptors',
    fingerprint: 'None',
    featureSelection: 'BorutaPy',
    tuning: 'Optuna',
    threshold: '0.5',
    homeDetails: [
      'Gradient Boosting classifier',
      '6 selected molecular descriptors',
      '189-compound in vivo corrosion dataset',
      'Test AUC 0.80, accuracy 0.79',
    ],
    metrics: [
      { label: 'Accuracy', train: '0.77', test: '0.79' },
      { label: 'AUC', train: '0.84', test: '0.80' },
      { label: 'F1', train: '0.73', test: '0.75' },
      { label: 'Recall', train: '0.72', test: '0.75' },
      { label: 'Precision', train: '0.73', test: '0.75' },
    ],
    hyperparameters: [
      'n_estimators: 50',
      'learning_rate: 0.01',
      'max_depth: 1',
      'subsample: 0.7',
      'min_samples_leaf: 10',
      'max_features: sqrt',
    ],
    documentationHref: `${docsBase}model-gb-corrosion-in-vivo`,
  },
  'GB Irritation (in vitro)': {
    modelName: 'GB Irritation (in vitro)',
    slug: 'model-gb-irritation-in-vitro',
    name: 'Skin Irritation 3D - In Vitro',
    badge: 'In Vitro',
    endpoint: 'Irritation',
    assayType: 'In vitro 3D',
    algorithm: 'Gradient Boosting classifier',
    dataset: 'In vitro 3D irritation dataset',
    observations: 208,
    classBalance: ['Positive (1): 116 / 0.56', 'Negative (0): 92 / 0.44'],
    descriptors: '7 molecular descriptors',
    fingerprint: 'None',
    featureSelection: 'RFE',
    tuning: 'Optuna',
    threshold: '0.5',
    homeDetails: [
      'Gradient Boosting classifier',
      '7 descriptors including HOMO and HL-gap',
      '208-compound in vitro 3D irritation dataset',
      'Test AUC 0.81, accuracy 0.79',
    ],
    metrics: [
      { label: 'Accuracy', train: '0.79', test: '0.79' },
      { label: 'AUC', train: '0.91', test: '0.81' },
      { label: 'F1', train: '0.78', test: '0.79' },
      { label: 'Recall', train: '0.79', test: '0.79' },
      { label: 'Precision', train: '0.81', test: '0.81' },
    ],
    hyperparameters: [
      'n_estimators: 55',
      'learning_rate: 0.017030189588951777',
      'max_depth: 1',
      'subsample: 0.5755346887476093',
      'min_samples_leaf: 5',
      'max_features: sqrt',
    ],
    documentationHref: `${docsBase}model-gb-irritation-in-vitro`,
  },
  'Ensemble Irritation (rabbit, in vivo)': {
    modelName: 'Ensemble Irritation (rabbit, in vivo)',
    slug: 'model-ensemble-irritation-rabbit',
    name: 'Skin Irritation - In Vivo Rabbit',
    badge: 'In Vivo',
    endpoint: 'Irritation',
    assayType: 'In vivo rabbit',
    algorithm: 'Stacking ensemble',
    dataset: 'In vivo rabbit irritation dataset',
    observations: 857,
    classBalance: ['Positive and negative cases from the 857-compound rabbit irritation dataset'],
    descriptors: '21 registered runtime descriptors',
    fingerprint: 'None',
    featureSelection: 'Individual model feature sets',
    tuning: 'Individual base model tuning + ensemble combination',
    threshold: '0.5',
    homeDetails: [
      'Stacking ensemble: XGB, RF, SVM, DT, KNN to LR',
      '21 registered runtime descriptors',
      '857-compound in vivo rabbit irritation dataset',
      'Stacking AUC 0.89, accuracy 0.8488',
    ],
    metrics: [
      { label: 'Accuracy', stacking: '0.8488', voting: '0.8449' },
      { label: 'AUC', stacking: '0.89', voting: '0.91' },
      { label: 'F1', stacking: '0.819', voting: '0.818' },
      { label: 'Recall', stacking: '0.819', voting: '0.833' },
      { label: 'Precision', stacking: '0.819', voting: '0.803' },
    ],
    hyperparameters: [
      'XGBoost: n_estimators 87, max_depth 4, learning_rate 0.075',
      'Random Forest: n_estimators 168, max_depth 5',
      'Decision Tree: max_depth 5',
      'KNN: n_neighbors 5',
      'SVM: C 0.366, kernel rbf',
    ],
    note: 'The final document lists 22 descriptors, but the runtime model and selected feature file use 21 descriptors. The runtime feature list is treated as authoritative.',
    documentationHref: `${docsBase}model-ensemble-irritation-rabbit`,
  },
  'XGB Phototox Chemico': {
    modelName: 'XGB Phototox Chemico',
    slug: 'model-xgb-phototox-chemico',
    name: 'Phototoxicity - In Chemico',
    badge: 'In Chemico',
    endpoint: 'Phototoxicity',
    assayType: 'In chemico',
    algorithm: 'XGBoost classifier',
    dataset: 'In chemico phototoxicity dataset',
    observations: 162,
    classBalance: ['Positive (1): 87 / 0.54', 'Negative (0): 75 / 0.46'],
    descriptors: '48 molecular descriptors',
    fingerprint: 'MACCS (166 bits)',
    featureSelection: 'Lasso',
    tuning: 'Optuna + Manual Fine-Tuning',
    threshold: '0.45',
    homeDetails: [
      'XGBoost classifier',
      '48 descriptors + MACCS fingerprint',
      '162-compound in chemico phototoxicity dataset',
      'Test AUC 0.87, accuracy 0.85',
    ],
    metrics: [
      { label: 'Accuracy', train: '0.89', test: '0.85' },
      { label: 'AUC', train: '0.90', test: '0.87' },
      { label: 'F1', train: '-', test: '0.85' },
      { label: 'Recall', train: '-', test: '0.78' },
      { label: 'Precision', train: '-', test: '0.93' },
    ],
    hyperparameters: [
      'n_estimators: 43',
      'max_depth: 2',
      'learning_rate: 0.16578669494542625',
      'subsample: 0.5092836919084545',
      'colsample_bytree: 0.5042728271034018',
      'min_child_weight: 2.053351654798258',
      'gamma: 0.013166513512511302',
      'reg_lambda: 1.9936713966098674',
      'reg_alpha: 0.010039348654324752',
    ],
    documentationHref: `${docsBase}model-xgb-phototox-chemico`,
  },
  'XGB Phototox 3T3': {
    modelName: 'XGB Phototox 3T3',
    slug: 'model-xgb-phototox-3t3',
    name: 'Phototoxicity 3T3 - In Vitro',
    badge: 'In Vitro',
    endpoint: 'Phototoxicity',
    assayType: 'In vitro 3T3',
    algorithm: 'XGBoost classifier',
    dataset: 'In vitro 3T3 NRU phototoxicity dataset',
    observations: 396,
    classBalance: ['Positive (1): 180 / 0.45', 'Negative (0): 216 / 0.55'],
    descriptors: '52 molecular descriptors',
    fingerprint: 'Atom Pair (512 bits)',
    featureSelection: 'XGB Importance',
    tuning: 'Randomized Search + Manual Fine-Tuning',
    threshold: '0.5',
    homeDetails: [
      'XGBoost classifier',
      '52 descriptors + AtomPair fingerprint',
      '396-compound in vitro 3T3 NRU dataset',
      'Test AUC 0.86, accuracy 0.83',
    ],
    metrics: [
      { label: 'Accuracy', train: '0.87', test: '0.83' },
      { label: 'AUC', train: '0.90', test: '0.86' },
      { label: 'F1', train: '-', test: '0.81' },
      { label: 'Recall', train: '-', test: '0.83' },
      { label: 'Precision', train: '-', test: '0.79' },
    ],
    hyperparameters: [
      'n_estimators: 70',
      'max_depth: 4',
      'learning_rate: 0.15',
      'subsample: 0.8',
      'colsample_bytree: 0.7',
      'min_child_weight: 1',
      'gamma: 0',
      'reg_lambda: 1',
      'reg_alpha: 0',
    ],
    documentationHref: `${docsBase}model-xgb-phototox-3t3`,
  },
  'XGB Phototox 3D': {
    modelName: 'XGB Phototox 3D',
    slug: 'model-xgb-phototox-3d',
    name: 'Phototoxicity 3D - In Vitro',
    badge: 'In Vitro',
    endpoint: 'Phototoxicity',
    assayType: 'In vitro 3D',
    algorithm: 'XGBoost classifier',
    dataset: 'In vitro 3D reconstructed-tissue phototoxicity dataset',
    observations: 101,
    classBalance: ['Positive (1): 45 / 0.45', 'Negative (0): 56 / 0.55'],
    descriptors: '24 molecular descriptors',
    fingerprint: 'RDKit (512 bits)',
    featureSelection: 'Lasso',
    tuning: 'Optuna + Manual Fine-Tuning',
    threshold: '0.5',
    homeDetails: [
      'XGBoost classifier',
      '24 descriptors + RDKit fingerprint',
      '101-compound in vitro 3D phototoxicity dataset',
      'Test AUC 0.89, accuracy 0.85',
    ],
    metrics: [
      { label: 'Accuracy', train: '0.86', test: '0.85' },
      { label: 'AUC', train: '0.90', test: '0.89' },
      { label: 'F1', train: '-', test: '0.84' },
      { label: 'Recall', train: '-', test: '0.88' },
      { label: 'Precision', train: '-', test: '0.80' },
    ],
    hyperparameters: [
      'n_estimators: 40',
      'max_depth: 3',
      'learning_rate: 0.23173378419139737',
      'subsample: 0.61309321921697',
      'colsample_bytree: 0.6438013267880465',
      'min_child_weight: 1.6133996725142632',
      'gamma: 0.7827579345694584',
      'reg_lambda: 0.24087072476030794',
      'reg_alpha: 0.013282230490271123',
    ],
    documentationHref: `${docsBase}model-xgb-phototox-3d`,
  },
  'SVM Phototox (in vivo)': {
    modelName: 'SVM Phototox (in vivo)',
    slug: 'model-svm-phototox-in-vivo',
    name: 'Phototoxicity - In Vivo',
    badge: 'In Vivo',
    endpoint: 'Phototoxicity',
    assayType: 'In vivo',
    algorithm: 'SVM pipeline',
    dataset: 'In vivo phototoxicity dataset',
    observations: 35,
    classBalance: ['Positive (1): 20 / 57%', 'Negative (0): 15 / 43%'],
    descriptors: 'RDKit descriptors + HOMO/LUMO candidates',
    fingerprint: 'AtomPair_128 candidate fingerprints',
    featureSelection: 'L1 Logistic Regression, max 30 selected features',
    tuning: 'Optuna 30 trials + RandomizedSearch',
    threshold: '0 (default)',
    homeDetails: [
      'SVM pipeline with embedded feature selection',
      '223 inputs with 30 selected internally',
      '35-compound in vivo phototoxicity dataset',
      'Test AUC 0.75, accuracy 0.8571',
    ],
    metrics: [
      { label: 'Accuracy', train: '1.00', test: '0.8571' },
      { label: 'AUC', train: '1.00', test: '0.75' },
      { label: 'F1', train: '1.00', test: '0.80' },
      { label: 'Recall', train: '1.00', test: '0.75' },
      { label: 'Precision', train: '1.00', test: '1.00' },
    ],
    hyperparameters: [
      'selector__estimator__C: 0.6931745808791483',
      'selector__max_features: 30',
      'svc__kernel: linear',
      'svc__C: 0.8803913550242988',
      'svc__class_weight: None',
    ],
    note: 'The validation set is small, so one classification change has a large effect on reported scores.',
    documentationHref: `${docsBase}model-svm-phototox-in-vivo`,
  },
};

export const modelDocumentationList = Object.values(modelDocumentationByName);

export function getModelDocumentation(modelName: string): ModelDocumentation | undefined {
  return modelDocumentationByName[modelName];
}
