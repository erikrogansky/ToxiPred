"""
Script to create database tables.
Run this to initialize or recreate the database schema.
"""
from app.db import engine, Base
from app.models.job_result import JobResult

def create_tables():
    """Drop all tables and recreate them with the new schema."""
    print("Dropping existing tables...")
    Base.metadata.drop_all(bind=engine)
    
    print("Creating tables with new schema...")
    Base.metadata.create_all(bind=engine)
    
    print("Database tables created successfully!")

if __name__ == "__main__":
    create_tables()
