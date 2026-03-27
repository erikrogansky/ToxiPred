#!/usr/bin/env python3
"""
Demo prediction management script.

Usage:
  # List current demos
  python manage_demos.py list

  # Seed demos (runs predictions and marks them as demos)
  python manage_demos.py seed

  # Clear all demo predictions
  python manage_demos.py clear

  # Replace all demos (clear + seed)
  python manage_demos.py replace

Environment:
  DATABASE_URL  – PostgreSQL or SQLite connection string
                  (default: sqlite:///./toxipred.db)
"""
import sys
import os

# Allow running from repo root or from server-toxipred/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from datetime import datetime, timedelta
from app.db import SessionLocal
from app.models.job_result import JobResult


# ────────────────────────────────────────────────────
# DEMO DEFINITIONS
# Edit this list to change what appears on the Demos page.
# Each entry needs:
#   - query: SMILES, CAS number, or trivial name
#   - model: model name exactly as registered (e.g. "XGBoost", "Ensamble")
#   - demo_title: short title shown on the card
#   - demo_description: 1-2 sentence description for the card
# ────────────────────────────────────────────────────
DEMO_COMPOUNDS = [
    {
        "query": "Aspirin",
        "model": "XGBoost",
        "demo_title": "Aspirin — In Vitro Phototoxicity",
        "demo_description": "A common NSAID analysed for phototoxic potential using the XGBoost model. See how molecular descriptors contribute to the prediction.",
    },
    {
        "query": "Caffeine",
        "model": "Ensamble",
        "demo_title": "Caffeine — In Chemico Photo-irritation",
        "demo_description": "Caffeine assessed for photo-irritation potential with the Ensemble classifier. Explore the full feature importance breakdown.",
    },
    {
        "query": "Ketoprofen",
        "model": "XGBoost",
        "demo_title": "Ketoprofen — Known Phototoxicant",
        "demo_description": "A well-known phototoxic NSAID. This demo shows how ToxiPred identifies phototoxic compounds and which features drive the prediction.",
    },
]


def list_demos():
    db = SessionLocal()
    try:
        demos = db.query(JobResult).filter(JobResult.is_demo == True).all()
        if not demos:
            print("No demo predictions found.")
            return
        for d in demos:
            pred = d.payload.get("prediction") if d.payload else None
            print(f"  [{d.id}] {d.demo_title or d.trivial_name or d.name}")
            print(f"    Model: {d.model} | Formula: {d.formula} | Prediction: {pred}")
            print()
    finally:
        db.close()


def clear_demos():
    db = SessionLocal()
    try:
        count = db.query(JobResult).filter(JobResult.is_demo == True).delete()
        db.commit()
        print(f"Cleared {count} demo prediction(s).")
    finally:
        db.close()


def seed_demos():
    """Submit predictions through the Celery pipeline and mark them as demos."""
    from app.services.chem_resolver import resolve_chemical
    from app.workers.celery_app import celery_app

    db = SessionLocal()
    try:
        existing = db.query(JobResult).filter(JobResult.is_demo == True).count()
        if existing:
            print(f"Warning: {existing} demo(s) already exist. Use 'replace' to clear first.")
            return

        for compound in DEMO_COMPOUNDS:
            print(f"Processing: {compound['query']} ({compound['model']})...")

            chem = resolve_chemical(compound["query"])
            if chem is None:
                print(f"  ERROR: Could not resolve '{compound['query']}'. Skipping.")
                continue

            # Send the prediction task synchronously
            task = celery_app.send_task(
                "toxipred.predict",
                args=[compound["model"], chem.smiles],
                kwargs={
                    "display_name": chem.name,
                    "formula": chem.formula,
                    "input_query": chem.query,
                    "input_type": chem.source,
                    "trivial_name": chem.trivial_name,
                    "other_names": chem.other_names,
                },
            )

            print(f"  Task {task.id} submitted. Waiting for result...")
            result = task.get(timeout=120)

            if not result:
                print(f"  ERROR: Task returned no result. Skipping.")
                continue

            # Find the job result in DB (the worker saves it) or create it
            job = db.get(JobResult, task.id)
            if job is None:
                job = JobResult.from_payload(task.id, result, ttl_days=365 * 100)
                db.add(job)

            # Mark as demo
            job.is_demo = True
            job.demo_title = compound["demo_title"]
            job.demo_description = compound["demo_description"]
            # Set a far-future expiry so it never expires
            job.expires_at = datetime.utcnow() + timedelta(days=365 * 100)

            db.commit()
            pred = result.get("prediction")
            print(f"  Done: {compound['demo_title']} (prediction={pred})")

        print("\nDemo seeding complete.")
    finally:
        db.close()


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    cmd = sys.argv[1]
    if cmd == "list":
        list_demos()
    elif cmd == "seed":
        seed_demos()
    elif cmd == "clear":
        clear_demos()
    elif cmd == "replace":
        clear_demos()
        seed_demos()
    else:
        print(f"Unknown command: {cmd}")
        print(__doc__)
        sys.exit(1)


if __name__ == "__main__":
    main()
