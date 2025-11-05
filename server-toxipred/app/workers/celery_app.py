# app/workers/celery_app.py
import os
from celery import Celery

BROKER  = os.getenv("CELERY_BROKER_URL",  "redis://toxipred-redis:6379/0")
BACKEND = os.getenv("CELERY_RESULT_BACKEND", "redis://toxipred-redis:6379/0")

# include= will let Celery import the module when it’s ready
celery_app = Celery(
    "toxipred",
    broker=BROKER,
    backend=BACKEND,
    include=["app.workers.tasks"],
)

celery_app.conf.update(
    task_track_started=True,
    task_default_queue=os.getenv("CELERY_DEFAULT_QUEUE", "default"),
)

# IMPORTANT: remove any manual `import app.workers.tasks` from here.