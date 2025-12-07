from datetime import datetime, timedelta
import secrets
import hashlib
from sqlalchemy import Column, String, DateTime, Integer, ForeignKey
from sqlalchemy.orm import relationship
from app.db import Base


def generate_share_token() -> str:
    """Generate a secure random token for sharing."""
    return secrets.token_urlsafe(16)


def hash_password(password: str) -> str:
    """Hash password using SHA-256."""
    return hashlib.sha256(password.encode()).hexdigest()


def verify_password(password: str, hashed: str) -> bool:
    """Verify password against hash."""
    return hash_password(password) == hashed


def generate_password() -> str:
    """Generate a random 16-character alphanumeric password."""
    return secrets.token_urlsafe(12)  # produces 16 chars


class SharedJob(Base):
    __tablename__ = "shared_jobs"

    id = Column(Integer, primary_key=True, autoincrement=True)
    share_token = Column(String, unique=True, nullable=False, index=True)
    password_hash = Column(String, nullable=False)
    source_job_id = Column(String, ForeignKey("job_results.id"), nullable=False)
    reference_count = Column(Integer, nullable=False, default=1)
    
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)
    expires_at = Column(DateTime, nullable=False, index=True)

    # Relationship to JobResult
    source_job = relationship("JobResult", backref="shared_links")

    @staticmethod
    def create_share(job_id: str, custom_password: str | None = None, ttl_days: int = 30) -> tuple["SharedJob", str]:
        """
        Create a new shared job entry.
        If custom_password is provided, use it; otherwise generate one.
        Returns the SharedJob instance and the plain-text password.
        """
        password = custom_password if custom_password else generate_password()
        now = datetime.utcnow()
        
        shared = SharedJob(
            share_token=generate_share_token(),
            password_hash=hash_password(password),
            source_job_id=job_id,
            reference_count=1,
            created_at=now,
            expires_at=now + timedelta(days=ttl_days),
        )
        return shared, password

    def verify(self, password: str) -> bool:
        """Verify the password for this share."""
        return verify_password(password, self.password_hash)

    def increment_refs(self) -> None:
        """Increment reference count when someone imports this share."""
        self.reference_count += 1

    def decrement_refs(self) -> int:
        """
        Decrement reference count.
        Returns the new count.
        """
        self.reference_count = max(0, self.reference_count - 1)
        return self.reference_count
