FROM mambaorg/micromamba:1.5.8
SHELL ["/bin/bash", "-lc"]

# Create env from conda-forge (pin NumPy 1.26 to match RDKit wheels)
RUN micromamba create -y -n toxipred -c conda-forge \
    python=3.11 \
    numpy=1.26 \
    scikit-learn=1.5.0 \
    rdkit \
    shap=0.45.1 \
    numba=0.59 \
    scipy=1.11 \
    pip \
 && micromamba clean -a -y

WORKDIR /app
COPY app ./app
COPY pyproject.toml ./

# Install the pip-only bits inside the env
# (leave ML/SCI libs to conda to avoid ABI mismatches)
RUN micromamba run -n toxipred python -m pip install --no-cache-dir \
  "fastapi>=0.115" \
  "uvicorn[standard]>=0.30" \
  "celery>=5.4" \
  "redis>=5.0" \
  "joblib>=1.3" \
  "xgboost>=2.0" \
  "dill>=0.3.8" \
  "pydantic>=2.7" \
  "requests>=2.31.0" \
  "flower>=2.0" \
  "SQLAlchemy>=2.0" \
  "psycopg2-binary>=2.9"

EXPOSE 8000
CMD ["micromamba","run","-n","toxipred","uvicorn","app.main:app","--host","0.0.0.0","--port","8000","--reload"]
