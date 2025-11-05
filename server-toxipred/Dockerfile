# server-toxipred/Dockerfile
FROM mambaorg/micromamba:1.5.8
SHELL ["/bin/bash", "-lc"]

# Create env from conda-forge (pins avoid the NumPy 2.x/RDKit ABI drama)
RUN micromamba create -y -n toxipred -c conda-forge \
    python=3.11 \
    numpy=1.26 \
    scikit-learn=1.4 \
    rdkit \
    pip \
 && micromamba clean -a -y

WORKDIR /app
COPY app ./app
COPY pyproject.toml ./

# Install the pip-only bits inside the env
RUN micromamba run -n toxipred python -m pip install --no-cache-dir \
  "fastapi>=0.115" \
  "uvicorn[standard]>=0.30" \
  "celery>=5.4" \
  "redis>=5.0" \
  "joblib>=1.3" \
  "xgboost>=2.0" \
  "dill>=0.3.8" \
  "pydantic>=2.7" \
  "flower>=2.0"

EXPOSE 8000
CMD ["micromamba","run","-n","toxipred","uvicorn","app.main:app","--host","0.0.0.0","--port","8000","--reload"]
