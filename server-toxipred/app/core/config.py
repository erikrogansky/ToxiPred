from pathlib import Path
import os

BASE_DIR = Path(__file__).resolve().parents[1]
MODELS_DIR = Path(os.getenv("TOXIPRED_MODELS_DIR", BASE_DIR / "models"))

REDIS_URL      = os.getenv("REDIS_URL", "redis://redis:6379/0")
CELERY_BROKER  = os.getenv("CELERY_BROKER", REDIS_URL)
CELERY_BACKEND = os.getenv("CELERY_BACKEND", "redis://redis:6379/1")
PROGRESS_DB    = int(os.getenv("PROGRESS_DB", "2"))
