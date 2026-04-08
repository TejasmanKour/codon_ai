import numpy as np
from sklearn.linear_model import LinearRegression
import pickle

# Dummy training data (you can replace later with real data)
X = np.array([
    [300, 0.45, 100],
    [450, 0.55, 150],
    [600, 0.60, 200],
    [250, 0.40, 80]
])

y = np.array([0.7, 0.85, 0.92, 0.6])  # expression scores

model = LinearRegression()
model.fit(X, y)

# Save model
with open("ai_engine/ml/model.pkl", "wb") as f:
    pickle.dump(model, f)

print("Model trained and saved!")