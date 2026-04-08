import os
import pickle
from ai_engine.ml.features import extract_features

def load_model():
    base_dir = os.path.dirname(__file__)
    model_path = os.path.join(base_dir, "model.pkl")

    with open(model_path, "rb") as f:
        return pickle.load(f)

# Load model once
model = load_model()

def predict_expression(sequence):
    features = extract_features(sequence)
    score = model.predict(features)[0]
    return round(float(score), 3)