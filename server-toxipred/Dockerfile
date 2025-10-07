# backend/Dockerfile
FROM python:3.11-slim

# 1) Fast Python installer
RUN pip install --no-cache-dir uv

WORKDIR /app
COPY pyproject.toml ./
# 2) Produce a pinned requirements file (cached layer)
RUN uv pip compile pyproject.toml -o requirements.txt
# 3) Install from the compiled lock
RUN uv pip install --system -r requirements.txt

# 4) Your code
COPY app ./app
EXPOSE 8000
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", "--reload"]
