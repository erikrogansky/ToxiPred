import pubchempy as pcp
from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors
import requests
import streamlit as st
import joblib
import pandas as pd
import json
import xgboost
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from streamlit import columns

features_list = ['MaxAbsEStateIndex', 'MaxEStateIndex', 'MinAbsEStateIndex', 'MinEStateIndex', 'qed', 'SPS', 'MolWt', 'HeavyAtomMolWt', 'ExactMolWt', 'NumValenceElectrons', 'NumRadicalElectrons', 'MaxPartialCharge', 'MinPartialCharge', 'MaxAbsPartialCharge', 'MinAbsPartialCharge', 'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3', 'BCUT2D_MWHI', 'BCUT2D_MWLOW', 'BCUT2D_CHGHI', 'BCUT2D_CHGLO', 'BCUT2D_LOGPHI', 'BCUT2D_LOGPLOW', 'BCUT2D_MRHI', 'BCUT2D_MRLOW', 'AvgIpc', 'BalabanJ', 'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v', 'Chi4n', 'Chi4v', 'HallKierAlpha', 'Ipc', 'Kappa1', 'Kappa2', 'Kappa3', 'LabuteASA', 'PEOE_VSA1', 'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'TPSA', 'EState_VSA1', 'EState_VSA10', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'VSA_EState1', 'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 'VSA_EState4', 'VSA_EState5', 'VSA_EState6', 'VSA_EState7', 'VSA_EState8', 'VSA_EState9', 'FractionCSP3', 'HeavyAtomCount', 'NHOHCount', 'NOCount', 'NumAliphaticCarbocycles', 'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAmideBonds', 'NumAromaticCarbocycles', 'NumAromaticHeterocycles', 'NumAromaticRings', 'NumAtomStereoCenters', 'NumBridgeheadAtoms', 'NumHAcceptors', 'NumHDonors', 'NumHeteroatoms', 'NumHeterocycles', 'NumRotatableBonds', 'NumSaturatedCarbocycles', 'NumSaturatedHeterocycles', 'NumSaturatedRings', 'NumUnspecifiedAtomStereoCenters', 'Phi', 'RingCount', 'MolLogP', 'MolMR', 'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 'fr_Ar_COO', 'fr_Ar_N', 'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO', 'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_Ndealkylation1', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde', 'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide', 'fr_aniline', 'fr_aryl_methyl', 'fr_benzene', 'fr_bicyclic', 'fr_ester', 'fr_ether', 'fr_guanido', 'fr_halogen', 'fr_imidazole', 'fr_ketone', 'fr_ketone_Topliss', 'fr_lactone', 'fr_methoxy', 'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 'fr_nitro_arom_nonortho', 'fr_para_hydroxylation', 'fr_phenol', 'fr_phenol_noOrthoHbond', 'fr_phos_acid', 'fr_phos_ester', 'fr_piperzine', 'fr_priamide', 'fr_pyridine', 'fr_quatN', 'fr_sulfide', 'fr_thiocyan', 'fr_thiophene', 'fr_unbrch_alkane']


def get_smiles_from_cas(cas_number):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{cas_number}/smiles"
    response = requests.get(url)

    if response.status_code == 200:
        smiles = response.text.strip()
        return smiles
    else:
        raise ValueError(f"Could not fetch SMILES for CAS number {cas_number}. HTTP Status Code: {response.status_code}")


def calculate_descriptors(smiles):
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None

        descriptor_names = [desc_name[0] for desc_name in Descriptors._descList]
        calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
        descriptors = calculator.CalcDescriptors(mol)
        descriptor_dict = dict(zip(descriptor_names, descriptors))

        return descriptor_dict
    except Exception:
        return None


def get_name_from_cas(cas_number):
    try:
        compounds = pcp.get_compounds(cas_number, 'name')
        if compounds and compounds[0].iupac_name:
            return compounds[0].iupac_name
        elif compounds and compounds[0].synonyms:
            return compounds[0].synonyms[0]  # fallback
        else:
            return "Substance name not found"
    except Exception:
        return "Error retrieving substance name"


model_files = {
    "XGBoost": "best_xgb_model.pkl",
    "Random Forest": "best_rf_model.pkl",
    "Logistic Regression": "best_logreg_model.pkl",
    "SVC": "best_svc_model.pkl",
    "Gradient Boosting": "best_gb_model.pkl"
}

scaler_files = {
    model: file.replace(".pkl", "_scaler.pkl")
    for model, file in model_files.items()
}

scalers = {}

feature_files = {
    model: file.replace(".pkl", "_features.json")
    for model, file in model_files.items()
}

models = {}
selected_features_dict = {}

for model_name, model_file in model_files.items():
    try:
        models[model_name] = joblib.load('models/' + model_file)
        with open('models/' + feature_files[model_name], "r") as f:
            selected_features_dict[model_name] = json.load(f)

        scaler_path = 'models/' + scaler_files[model_name]
        try:
            scalers[model_name] = joblib.load(scaler_path)
        except FileNotFoundError:
            scalers[model_name] = None
    except Exception as e:
        st.warning(f"Failed to load model or features for {model_name}: {e}")

st.title("Irritation Prediction from CAS Number")

# User input
cas_number = st.text_input("Enter CAS Number:")

# Model selection
model_choice = st.selectbox("Select Model:", list(models.keys()))

if st.button("Predict"):
    if not cas_number:
        st.error("Please enter a CAS number.")
    else:
        try:
            smiles = get_smiles_from_cas(cas_number)
            if smiles:
                st.success(f"SMILES code: {smiles}")
                substance_name = get_name_from_cas(cas_number)
                st.write(f"Substance name: {substance_name}")

                # Calculate descriptors from SMILES
                descriptors = calculate_descriptors(smiles)
                if descriptors is None or len(descriptors) == 0:
                    st.error("Failed to compute molecular descriptors. Please check the CAS number.")
                    st.stop()

                descriptors_df = pd.DataFrame([descriptors])
                descriptors_df = descriptors_df[features_list]

                scaler = scalers.get(model_choice, None)
                if scaler is not None:
                    descriptors_df = pd.DataFrame(
                        scaler.transform(descriptors_df),
                        columns=descriptors_df.columns
                    )

                selected_features = selected_features_dict[model_choice]

                descriptors_df = descriptors_df[selected_features]
                descriptors_df = descriptors_df.fillna(0)

                model = models[model_choice]

                prediction_proba = model.predict_proba(descriptors_df)[:, 1]

                prediction = model.predict(descriptors_df)[0]

                if prediction_proba[0] is not None:
                    st.write(f"Prediction Probability of Irritant: {prediction_proba[0]:.2f}")
                else:
                    st.write("Probability not available for this model.")

                st.write(f"Predicted Class: {'Irritant' if prediction == 1 else 'Non-Irritant'}")
            else:
                st.error("Failed to retrieve SMILES for the given CAS number.")
        except Exception as e:
            st.error(f"An error occurred: {e}")
